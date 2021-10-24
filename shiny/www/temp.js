let db;
Shiny.addCustomMessageHandler('set-user', function(msg) {
	$('#authentication-user').text(msg.user);
	db = firebase.firestore();
});

function logout(){
	Shiny.setInputValue('load-auth-firebaseLogout', 1, {priority: 'event'});
}

function upgrade(){
	const docRef = await db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.collection('checkout_sessions')
		.add({
			price: 'price_1Jo2cULGmSWfyZoW6RUicoX6',
			success_url: window.location.origin,
			cancel_url: window.location.origin,
		});

	// Wait for the CheckoutSession to get attached by the extension
	docRef.onSnapshot((snap) => {
		const { error, url } = snap.data();
		if (error) {
			// Show an error to your customer and 
			// inspect your Cloud Function logs in the Firebase console.
			alert(`An error occured: ${error.message}`);
		}
		if (url) {
			// We have a Stripe Checkout URL, let's redirect.
			window.location.assign(url);
		}
	});	
}
