##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#--------------------------------------------------------------------------------
# This JS code is responsible for activating and desactivating the class `active`,
# defined on `scss/components/_dropdownmenu.scss` when clicking a dropdown button.
# It also has logic implemented so that when clicking a different dropdown item,
# the active one gets desactivated.

## second version of Carson (see email 9jul23)
dropdown.jsCode.SAVE <- function(id) {
  paste0("
$('#", id, "').on('click', function(){
  if($('#", id, "').hasClass('active')){
    $('#", id, "').toggleClass('active');
    return 0;
  }
  $('.dropdown-button.active').toggleClass('active');
  if(!$('#", id, "').hasClass('active')){
    $('#", id, "').toggleClass('active');
  }
})
function restoreDropdownMenu() {
  var menuBody = $('body > .dropdown-menu-body');
  if (menuBody.length > 0) {
    menuBody.removeClass('dropdown-menu-body');
    var id = menuBody.attr('data-toggle-id');
    $('#' + id).parent().append(menuBody.detach());
  }
}
$('#", id, "').on('shown.bs.dropdown', function () {
  restoreDropdownMenu();
  if (!this.id) {
    console.error('Expected dropdown toggle to have an id attribute');
  }
  setTimeout(function() {
    var menu = $(this).parent().find('.dropdown-menu');
    menu.addClass('dropdown-menu-body');
    menu.attr('data-toggle-id', this.id);
    $('body').append(menu.detach());
  },1);
});
$('#", id, "').on('hidden.bs.dropdown', function () {
  restoreDropdownMenu();
});
") ## end of paste0
}

## first version of Carson
dropdown.jsCode <- function(id) {
  paste0("
$('#", id, "').on('click', function(){
  if($('#", id, "').hasClass('active')){
    $('#", id, "').toggleClass('active');
    return 0;
  }
  $('.dropdown-button.active').toggleClass('active');
  if(!$('#", id, "').hasClass('active')){
    $('#", id, "').toggleClass('active');
  }
})
function restoreDropdownMenu() {
  var menuBody = $('body > .dropdown-menu-body');
  if (menuBody.length > 0) {
    menuBody.removeClass('dropdown-menu-body');
    var id = menuBody.attr('data-toggle-id');
    $('#' + id).parent().append(menuBody.detach());
  }
}
$('#", id, "').on('shown.bs.dropdown', function () {
  restoreDropdownMenu();
  if (!this.id) {
    console.error('Expected dropdown toggle to have an id attribute');
  }
  var menu = $(this).parent().find('.dropdown-menu');
  menu.addClass('dropdown-menu-body');
  menu.attr('data-toggle-id', this.id);
  $('body').append(menu.detach());
});
$('#", id, "').on('hidden.bs.dropdown', function () {
  restoreDropdownMenu();
});
") ## end of paste0
}

DropdownMenu <- function(..., size = "default", status = "default", icon = NULL, width = "auto", margin = "10px") {
  id <- bigdash:::make_id()
  tags$div(
    tags$a(
      class = paste0(
        "btn btn-", status, " action-button dropdown-button btn-circle", if (size == "default") {
          " "
        } else {
          paste0("-", size, " ")
        },
        if (size != "default") paste0("btn-", size)
      ),
      id = id, type = "button", `data-bs-toggle` = "dropdown",
      `aria-expanded` = "false", list(icon)
    ),
    tags$div(
      class = "dropdown-menu",
      `aria-labelledby` = id,
      tryCatch(
        {
          lapply(list(...), function(x) {
            htmltools::tags$li(
              x,
              style = paste0(
                "width: ", htmltools::validateCssUnit(width),
                "; margin-left: ", htmltools::validateCssUnit(margin),
                "; margin-right: ", htmltools::validateCssUnit(margin), ";"
              )
            )
          })
        },
        error = function(w) {}
      )
    ),
    tags$script(HTML(dropdown.jsCode(id)))
  )
}

# same thing as DropdownMenu but without circle icon;
# currently used in the board.loading table
actionMenu <- function(..., size = "default", status = "default", icon = NULL, margin = "10px") {
  id <- bigdash:::make_id()
  tags$div(
    tags$a(
      class = paste0(
        "btn btn-outline-", status, " action-button dropdown-button", if (size == "default") {
          " "
        } else {
          paste0("-", size, " ")
        },
        if (size != "default") paste0("btn-", size)
      ),
      style = "border: none;",
      id = id, type = "button", `data-bs-toggle` = "dropdown",
      `aria-expanded` = "false", list(icon)
    ),
    tags$div(
      class = "dropdown-menu",
      `aria-labelledby` = id,
      tryCatch(
        {
          lapply(list(...), function(x) {
            htmltools::tags$li(
              x,
              style = paste0(
                "margin-left: ",
                htmltools::validateCssUnit(margin), "; margin-right: ",
                htmltools::validateCssUnit(margin), ";"
              )
            )
          })
        },
        error = function(w) {}
      )
    ),
    tags$script(HTML(dropdown.jsCode(id)))
  )
}
