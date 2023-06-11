const quit = () => {
    Shiny.setInputValue('quit', 1, {priority: 'event'});
};

const shinyproxy_logout = () => {
	window.location.assign("/logout");
};
