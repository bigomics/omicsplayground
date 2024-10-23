$(document).on('click', '[id$="toggle_button"]', function() {
      var $card = $(this).closest('.card');
      var isFullScreen = $card.attr('data-full-screen') === 'true';

      var elementId = $(this).attr('id');
      var prefix = elementId.split('toggle_button')[0];

      if (!isFullScreen) {
        // Entering full screen: add `bslib-full-screen-overlay` as sibling before card
        var $overlay = $('<div>', { id: 'bslib-full-screen-overlay' });
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
