var shouldWarnOnExit = false;

Shiny.addCustomMessageHandler('warnOnExit', function(message) {
    shouldWarnOnExit = message;
});

// Before the user closes the tab
window.addEventListener('beforeunload', function (e) {
if (shouldWarnOnExit) {
        var confirmationMessage = 'Dataset computation is running. Are you sure you want to leave? Leaving will cancel the computation.';
        (e || window.event).returnValue = confirmationMessage; // For older browsers
        return confirmationMessage; // For modern browsers
    }
});