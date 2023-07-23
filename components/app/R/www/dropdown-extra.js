// Carson solution so that dropdowns appear to be on top of the
// bslib cards
$('.dropdown-button').on('click', function(){
    var id = $(this).attr('id');
    if($('#' + id).hasClass('active')){
        $('#' + id).toggleClass('active');
        return 0;
    }
    $('.dropdown-button.active').toggleClass('active');
    if(!$('#' + id).hasClass('active')){
        $('#' + id).toggleClass('active');
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
  
  $('.dropdown-button').on('shown.bs.dropdown', function () {
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
  
  $('.dropdown-button').on('hidden.bs.dropdown', function () {
    restoreDropdownMenu();
  });