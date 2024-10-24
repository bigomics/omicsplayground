$(document).on('click', '[id$="toggle_button"]', function() {
      var $card = $(this).closest('.card');
      if ($card.length === 0) {
          $card = $(this).parent().next();
      }
      var isFullScreen = $card.attr('data-full-screen') === 'true';

      var elementId = $(this).attr('id');
      var prefix = elementId.split('toggle_button')[0];

      if (!isFullScreen) {
        // Entering full screen: add `bslib-full-screen-overlay` as sibling before card
        var $overlay = $('<div>', { id: 'bslib-full-screen-overlay' }).append(
          $('<a>', {
            id: prefix + 'toggle_button',
            class: 'bslib-full-screen-exit',
            tabindex: 0,
            html: 'Close <svg width=\"20\" height=\"20\" fill=\"currentColor\" class=\"bi bi-x-lg\" viewBox=\"0 0 16 16\"><path d=\"M2.146 2.854a.5.5 0 1 1 .708-.708L8 7.293l5.146-5.147a.5.5 0 0 1 .708.708L8.707 8l5.147 5.146a.5.5 0 0 1-.708.708L8 8.707l-5.146 5.147a.5.5 0 0 1-.708-.708L7.293 8 2.146 2.854Z\"></path></svg>'
          })
        );
        $card.before($overlay);  // Insert the overlay before the card
        $card.attr('data-full-screen', 'true');
        Shiny.setInputValue(prefix + 'fullscreen', true);
      } else {
        // Exiting full screen: remove the `bslib-full-screen-overlay`
        $card.prev('#bslib-full-screen-overlay').remove();  // Remove the overlay
        $card.attr('data-full-screen', 'false');
        Shiny.setInputValue(prefix + 'fullscreen', false);
      }
});
