##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# The four OPG AI menus, in settings-UI order. Every provider is described
# against this fixed vocabulary so the rest of the code can iterate menus
# without hard-coding names per provider.
.opg_ai_menus <- c("reports", "images", "copilot_deep", "copilot_balanced")

# Coerce a policy field (which may be NULL, a scalar, or a nested list from
# JSON) into a clean character vector with NAs and blanks dropped, so every
# downstream step can treat all fields uniformly.
.opg_ai_chr <- function(x) {
  if (is.null(x)) return(character(0))
  x <- unlist(x, use.names = FALSE)
  x[!is.na(x) & nzchar(x)]
}

# Load the OPG AI selection policy from JSON. A missing file yields an empty
# policy (menus then fall back to whatever the catalog offers); a wrong "kind"
# warns rather than stops, so a provider inventory supplied by mistake surfaces
# early instead of silently mis-filtering menus.
.opg_ai_read_policy <- function(path) {
  if (!file.exists(path)) {
    message("[GLOBAL] ai_model_policy.json not found, using empty AI selection policy")
    return(list())
  }
  policy <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(policy$kind, "opg_ai_selection_policy")) {
    warning(
      "[GLOBAL] ai_model_policy.json is expected to be an OPG AI selection policy, ",
      "not a provider inventory.",
      call. = FALSE
    )
  }
  policy
}

# Narrow and order model ids by one menu's policy: keep only `include` (when
# given), drop `exclude`, then float `prefer` entries to the front so the first
# surviving id becomes the menu's default selection.
.opg_ai_apply_menu_policy <- function(ids, menu_policy) {
  ids <- .opg_ai_chr(ids)
  if (length(ids) == 0L) return(ids)

  include <- .opg_ai_chr(menu_policy$include)
  if (length(include)) ids <- ids[ids %in% include]

  exclude <- .opg_ai_chr(menu_policy$exclude)
  if (length(exclude)) ids <- setdiff(ids, exclude)

  prefer <- .opg_ai_chr(menu_policy$prefer)
  if (length(prefer)) ids <- c(intersect(prefer, ids), setdiff(ids, prefer))

  unique(ids)
}

# Resolve a menu's models from omicsai's live catalog: take every known model
# for the provider that matches the menu's required capability, then apply the
# include/exclude/prefer policy. Keeps model facts owned by omicsai while OPG
# owns only the curation.
.opg_ai_catalog_menu <- function(provider, menu_policy) {
  capability <- .opg_ai_chr(menu_policy$capability)
  if (length(capability) == 0L) return(character(0))

  rows <- omicsai::ai_known_models(provider = provider, capability = capability)
  if (nrow(rows) == 0L) return(character(0))

  .opg_ai_apply_menu_policy(rows$id, menu_policy)
}

# Build all four menus for one provider. A menu carrying explicit `models`
# uses them verbatim; otherwise, when enabled, it is filled from the capability
# catalog. This lets a deployment pin exact models (e.g. bigomics) or defer to
# whatever omicsai currently offers (e.g. anthropic).
.opg_ai_provider_menus <- function(provider, policy) {
  provider_policy <- policy$providers[[provider]] %||% list()
  provider_menus <- provider_policy$menus %||% list()
  menu_policy <- policy$menus %||% list()

  stats::setNames(lapply(.opg_ai_menus, function(menu) {
    local_policy <- provider_menus[[menu]] %||% list()
    direct_models <- .opg_ai_chr(local_policy$models)

    if (length(direct_models)) {
      return(direct_models)
    }
    if (isTRUE(local_policy$enabled) || is.null(local_policy$enabled)) {
      return(.opg_ai_catalog_menu(
        provider = provider,
        menu_policy = utils::modifyList(menu_policy[[menu]] %||% list(),
                                        local_policy)
      ))
    }
    character(0)
  }), .opg_ai_menus)
}

# Build the provider -> menu -> models map for every enabled provider. This is
# the single structure the settings server reads to populate its per-provider
# model menus, computed once at startup.
.opg_ai_build_models <- function(policy, enabled_providers) {
  enabled_providers <- .opg_ai_chr(enabled_providers)
  models <- stats::setNames(vector("list", length(enabled_providers)),
                            enabled_providers)
  for (provider in enabled_providers) {
    models[[provider]] <- .opg_ai_provider_menus(provider, policy)
  }
  models
}

# Union of one menu's models across the given providers. Seeds the initial,
# provider-agnostic menu choices shown before the user picks a provider (the
# per-provider observer narrows them afterwards).
.opg_ai_menu_allowlist <- function(models, providers, menu) {
  entries <- unlist(lapply(providers, function(p) models[[p]][[menu]]),
                    use.names = FALSE)
  unique(entries[!is.na(entries) & nzchar(entries)])
}
