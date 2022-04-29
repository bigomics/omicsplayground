let db;
let pricing;
Shiny.addCustomMessageHandler('set-user', function(msg) {
	$('#authentication-user').text(msg.user);
	pricing = msg.pricing;
	if(msg.level == "premium"){
		// $('#authentication-upgrade').hide();  // really?
	}
});

Shiny.addCustomMessageHandler('manage-sub', function(msg) {
	window.location.assign(msg);
});

Shiny.addCustomMessageHandler('get-permissions', function(msg) {
	if(!db)
		db = firebase.firestore();
	
	db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.get()
		.then(function(doc){
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

Shiny.addCustomMessageHandler('get-subs', function(msg) {
	if(!db)
		db = firebase.firestore();

	db
		.collection('customers')
		.doc(firebase.auth().currentUser.uid)
		.collection('subscriptions')
		.where('status', 'in', ['active', 'inactive'])
		.onSnapshot(async (snapshot) => {
			const docs = [];

			snapshot.docs.forEach(function(doc){
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
			docs.forEach(function(doc) {
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

		}, function(error) {
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

function logout(){
    Shiny.setInputValue('auth-firebaseLogout', 1, {priority: 'event'});    
};

function quit(){
    Shiny.setInputValue('quit', 1, {priority: 'event'});  // trigger shiny input$quit
    // window.close();  // close window??
};

Shiny.addCustomMessageHandler('shinyproxy-logout', function(msg) {
    window.location.assign("/logout");
});


function show_plans(){
    Shiny.setInputValue('auth-firebaseUpgrade', 1, {priority: 'event'});    
};

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

function toggleEmail(){
	$('#emailLinkWrapper').toggle();
}

Shiny.addCustomMessageHandler('email-feedback', function(msg) {
	$('#emailFeedbackShow').html(msg.msg);

	setTimeout(function() {
		$('#emailFeedbackShow').html('');
	}, 5000);
});

function priceChange(name){
	if(name == 'monthly'){
		$('#yearlyCheck').prop('checked', false);
	} else {
		$('#monthlyCheck').prop('checked', false);
	}		
	if($('#yearlyCheck').prop('checked')){
		$('#starter-pricing').text('CHF49 / month');
		$('#premium-pricing').text('CHF490 / month');
	} else {
		$('#starter-pricing').text('CHF69 / month');
		$('#premium-pricing').text('CHF690 / month');
	}
}

function hideSub() {
    $('#logSub').hide();
    $('#logSubbed').show();
    $('#sever-reload-btn').show();
}

//function showSub() {
//	$('#logSub').hide();
//	$('#logSubbed').show();
//}

function sendLog() {
//	showSub();
        let msg  = $('#logMsg').val();
//	let user = $('#authentication-user').val();

	fetch(`log?msg=${encodeURIComponent(msg)}`)
		.then(res => {
			console.info(res);
			hideSub();
		})
	        .catch(error => {
			console.error(error);
			hideSub();
		})
}

function sendLog2(msg){
//	showSub();
//      let msg  = $('#logMsg').val();
//	let user = $('#authentication-user').val();

//    	fetch(`log?session=${id}&msg=${encodeURIComponent(msg)}`)
    	fetch(`log?msg=${encodeURIComponent(msg)}`)    
	        .then(res => {
		        console.info(res);
			hideSub();
		})
	        .catch(error => {
			console.error(error);
			hideSub();
		})
};

Shiny.addCustomMessageHandler('referral-input-error', function(msg) {
	$(`#${msg.target}`).addClass('error');
	$(`#${msg.target}`).after(`<small class='text-danger'>${msg.message}</small>`);

	setTimeout(() => {
		$(`#${msg.target}`).removeClass('error');
		$(`#${msg.target}`)
			.siblings('small')
			.remove();
	}, 5000);
});

Shiny.addCustomMessageHandler('referral-global-error', function(msg) {
	$('#referral-global-error').html(msg.message);
	setTimeout(() => {
		$('#referral-global-error').html('');
	}, 5000);
});

$(() => {
	$('.sidebar-content')
		.children()
		.each((index, el) => {
			if($(el).hasClass('collapse'))
				return;

			if(index == 0){
				$(el).show();
				return;
			} 
			
			$(el).hide();
		});
});

Shiny.addCustomMessageHandler('show-tabs', function(msg) {
	setTimeout(() => {
		$('.sidebar-content')
			.children()
			.each((index, el) => {
				if($(el).hasClass('collapse'))
					return;

				$(el).show();
			});
	}, 1000);
});
