// JS function to make the popovers
// dismissable when clcking outside the
// popover itself. 
// Called from ui-DropdownMenu.R
function makePopoverDismissible(id) {
  $(document).click(function (e) {
    var popover = $('#' + id);
    //popover.addClass("btn-active");
    if (!popover.is(e.target) && popover.has(e.target).length === 0 && $('.popover').has(e.target).length === 0) {
      popover.popover('hide');
      // popover.removeClass("btn-active");
    }
  });
}

// JS function to add and remove class `btn-active`
// when popover is shown or hidden
function addActionOnPopoverChange(popoverId) {
  $('#' + popoverId).on('hidden.bs.popover', function () {
    $(this).removeClass('btn-active');
  });
  $('#' + popoverId).on('shown.bs.popover', function () {
    $(this).addClass('btn-active');
  });
}
