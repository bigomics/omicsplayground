Shiny.addCustomMessageHandler('set-user', function(msg) {
	$('#authentication-user').text(msg.user);
});

function logout(){
	Shiny.setInputValue('load-auth-firebaseLogout', 1, {priority: 'event'});
}

function upgrade(){
	console.log("Coming soon")
}
