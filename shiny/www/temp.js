Shiny.addCustomMessageHandler('set-user', function(msg) {
	$('#authentication-user').text(msg.user);
});

function logout(){
    Shiny.setInputValue('load-auth-firebaseLogout', 1, {priority: 'event'});
};

function upgrade(){
    console.log("Coming soon");
    Shiny.setInputValue('load-auth-firebaseUpgrade', 1, {priority: 'event'});    
};

function sigstop(){
    Shiny.setInputValue('sigstop', 1, {priority: 'event'});  // end session
    window.close();  // close window??
};
