
#' Return available tools from omicsagentovi to be used as Ellmer
#' ToolDefs.
#' 
obi.get_tools <- function(pgx, group=NULL, names=NULL) {
  require(omicsagentovi)
  
  # 2. Mutable runtime env
  runtime <- new.env(parent = emptyenv())
  runtime$context <- omicsagentovi::context_set_pgx(
    omicsagentovi::RunContext(), pgx)
  runtime$run_state <- omicsagentovi::run_state_start(
    omicsagentovi::RunState())   # <-- "idle" → "running"
  
  # 3. Instantiate tools
  tools <- omicsagentovi::ovi_tools(
    session   = omicsagentovi::AgentSession(),
    run_state = omicsagentovi::RunState(),
    bindings  = omicsagentovi::RunBindings(),
    runtime   = runtime,
    group     = group,
    names     = names
  )

  return(tools)
}
