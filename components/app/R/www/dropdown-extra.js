// --------------------------------------------------------------------------------
//  This JS code is responsible for activating and desactivating the class `active`,
//  defined on `scss/components/_dropdownmenu.scss` when clicking a dropdown button.
//  It also has logic implemented so that when clicking a different dropdown item,
//  the active one gets desactivated.

// Carson solution so that dropdowns appear to be on top of the
// bslib cards
function initializeDropdown(dropdownId) {
  var dropdownSelector = '#' + dropdownId;

  $(dropdownSelector).on('click', function(){
      if($(dropdownSelector).hasClass('active')){
          $(dropdownSelector).toggleClass('active');
          return 0;
      }
      $('.dropdown-button.active').toggleClass('active');
      if(!$(dropdownSelector).hasClass('active')){
          $(dropdownSelector).toggleClass('active');
      }
  });
  
  $(dropdownSelector).on('shown.bs.dropdown', function () {
      restoreDropdownMenu();
      if (!this.id) {
          console.error('Expected dropdown toggle to have an id attribute');
      }
      var menu = $(this).parent().find('.dropdown-menu');
      menu.addClass('dropdown-menu-body');
      menu.attr('data-toggle-id', this.id);
      menu.detach();
      setTimeout(function() { $('body').append(menu); }, 0);
  });

  $(dropdownSelector).on('hidden.bs.dropdown', function () {
      restoreDropdownMenu();
  });

  function restoreDropdownMenu() {
      var menuBody = $('body > .dropdown-menu-body');
      if (menuBody.length > 0) {
          menuBody.removeClass('dropdown-menu-body');
          var id = menuBody.attr('data-toggle-id');
          $('#' + id).parent().append(menuBody.detach());
      }
  }
}
