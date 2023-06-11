let db;
let pricing;
let user;

Shiny.addCustomMessageHandler('get-permissions', (msg) => {
	if(!db)
		db = firebase.firestore();

	db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.get()
		.then((doc) => {
			let data = {
				href: window.location.href,
				id: doc.data().stripeId
			}
			Shiny.setInputValue(`${msg.ns}-stripeId`, data);
		});


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
		}, (error) => {
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

Shiny.addCustomMessageHandler('get-subs', (msg) => {
	if(!db)
		db = firebase.firestore();

	db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.collection('subscriptions')
		.where('status', 'in', ['active', 'inactive'])
		.onSnapshot(async (snapshot) => {
			const docs = [];

			snapshot.docs.forEach((doc) => {
				let docp = doc.data();
				docs.push({
					status: docp.status,
					cancel_at_period_end: docp.cancel_at_period_end,
					current_period_start: docp.current_period_start.seconds,
					current_period_end: docp.current_period_end.seconds,
				});
			})

			// no past subs
			if(docs.length == 0){
				$('#user-subs').html(
					'<span>No past subscriptions</span>'
				);
				return;
			}

			let content = '';
			docs.forEach((doc) => {
				let cl = 'warning';
				if(doc.status == 'active')
					cl = 'success';

				content += `<div class='thumbnail'>
					<div class='caption'>
						<div class='row'>
							<div class='col-md-3'>
								<br/>
								Monthly Subscription
							</div>
							<div class='col-md-3'>
								<span style='color:grey;'>Start</span><br/>
								${new Date(doc.current_period_start * 1000).toLocaleDateString()}
							</div>
							<div class='col-md-3'>
								<span style='color:grey;'>End</span><br/>
								${new Date(doc.current_period_end * 1000).toLocaleDateString()}
							</div>
							<div class='col-md-3'>
								<span style='color:grey;'>Status</span><br/>
								<span class='label label-${cl}'>${doc.status}</span>
							</div>
						</div>
					</div>
				</div>`;
			});

			$('#user-subs').html(content);

		}, (error) => {
			Shiny.setInputValue(
				msg.ns + '-subs',
				{
					success: false,
					response: error
				},
				{priority: 'event'}
			);
		});
});

async function upgrade_plan(){
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

const show_plans = () => {
    Shiny.setInputValue('auth-firebaseUpgrade', 1, {priority: 'event'});
  };