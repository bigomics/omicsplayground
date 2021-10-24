let db;
let pricing;
Shiny.addCustomMessageHandler('set-user', function(msg) {
	$('#authentication-user').text(msg.user);
	pricing = msg.pricing;
});

Shiny.addCustomMessageHandler('get-permissions', function(msg) {
	if(!db)
		db = firebase.firestore();

	db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.collection('subscriptions')
		.where('status', '==', 'active')
		.onSnapshot(async (snapshot) => {
			const doc = snapshot.docs[0];

			try {
				Shiny.setInputValue(
					msg.ns + '-permissions', 
					{ 
						success: true,
						response: doc.data() 
					},
					{priority: 'event'}
				);
			} catch (error) {
				Shiny.setInputValue(
					msg.ns + '-permissions', 
					{ 
						success: false,
						response: error
					},
					{priority: 'event'}
				);
			}
		}, function(error) {
			Shiny.setInputValue(
				msg.ns + '-permissions', 
				{ 
					success: false,
					response: error
				},
				{priority: 'event'}
			);
		});
});

function logout(){
	Shiny.setInputValue('load-auth-firebaseLogout', 1, {priority: 'event'});
}

async function upgrade(){
	const docRef = await db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.collection('checkout_sessions')
		.add({
			price: pricing,
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
