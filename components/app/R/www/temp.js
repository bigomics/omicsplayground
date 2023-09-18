let db;
let pricing;
let user;

//$(document).ready(function() {
$(document).on('shiny:connected', function() {
	
	// prevent browser back and forward buttons from working
	history.pushState(null, null, location.href);
	window.onpopstate = function () {
		history.go(1);
	};


    $(document).on('change', '.card-footer-checked', function(e) {
	// Set the "checked" property for all the card-footer-checked elements
	var isChecked = $(this).prop("checked");
        $(".card-footer-checked").prop("checked", isChecked);

	if ($(this).prop("checked") === true) {
	    $(".card-footer").show().animate({height: '3rem'}, 200);
	};

	if ($(this).prop("checked") === false) {
	    $(".card-footer").animate({height: '0px'}, 200, function() {
		$(this).hide();
	    });
	};

	$(window).resize();  // yikes...
	console.log('window.resize!');
    });

// call with: session$sendCustomMessage("window_resize", list(resize = TRUE))
    Shiny.addCustomMessageHandler('window_resize', function(message) {
        // console.log('hit', message.resize);
        if (message.resize) $(window).resize();
    })

    Shiny.addCustomMessageHandler('set-user', (msg) => {
//        $('#authentication-user').text(msg.user);
        user = msg.user;
	pricing = msg.pricing;
	if(msg.level == "premium"){
	    // $('#authentication-upgrade').hide();  // really?
	}
    });    

});  // end of on.shiny.connected


const unloadSidebar = () => {
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
        $('#sidebar-help-container').hide();
}

const sidebarClose = () => {
    if($('#sidebar-container').hasClass('sidebar-expanded')) {
	$('.sidebar-label').trigger('click');
    }
    $('#sidebar-help-container').hide();
}

const sidebarOpen = () => {
    if($('#sidebar-container').hasClass('sidebar-collapsed')) {
	$('.sidebar-label').trigger('click');
    }
    $('#sidebar-help-container').show();
}

const settingsClose = () => {
	if($('#settings-container').hasClass('settings-expanded'))
		$('#settings-container').trigger('mouseleave');

}

const settingsOpen = () => {
	if($('#settings-container').hasClass('settings-collapsed'))
		$('#settings-container').trigger('mouseenter');
}

const settingsLock = () => {
	if($('#settings-container').hasClass('settings-unlocked'))
		$('.settings-lock').trigger('click');
	if(!$('#settings-container').hasClass('settings-locked'))
		$('.settings-lock').trigger('click');
}

const settingsUnlock = () => {
	if($('#settings-container').hasClass('settings-locked'))
		$('.settings-lock').trigger('click');
}

$(() => {
	unloadSidebar();
});

$(function(){
        // init sequence: close sidebar, goto Welcome page and hide it's tab item
	setTimeout(() => {
		$('.sidebar-label').trigger('click');
		$('.sidebar-menu')
			.first()
			.trigger('click');

		$('.tab-sidebar')
			.first()
			.trigger('click');
			//.css('display', 'none');
		// on mouseover this does not work anymore, substitute by lock button option
	        //$('.settings-label').click()
	}, 250);

	const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
	const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));
})


Shiny.addCustomMessageHandler('manage-sub', (msg) => {
	window.location.assign(msg);
});

Shiny.addCustomMessageHandler('enableInfo', (data) => {
  Shiny.setInputValue(data.id, data.value);
});

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

const logout = () => {
	unloadSidebar();
	sidebarClose();
 	Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

const logoutInApp = () => {
	unloadSidebar();
        $(".tab-trigger[data-target='welcome-tab']").trigger('click');
	sidebarClose();
        Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

const quit = () => {
        Shiny.setInputValue('quit', 1, {priority: 'event'});
};

Shiny.addCustomMessageHandler('shinyproxy-logout', (msg) => {
        window.location.assign("/logout");
});

/* ********************* STRIPE HANDLERS **************************** */
const show_plans = () => {
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

Shiny.addCustomMessageHandler('email-feedback', function(msg) {
	$('#emailFeedbackShow').html(msg.msg);

	setTimeout(function() {
		$('#emailFeedbackShow').html('');
	}, 5000);
});

const priceChange = (name) => {
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


Shiny.addCustomMessageHandler('referral-input-error', (msg) => {
	$(`#${msg.target}`).addClass('error');
	$(`#${msg.target}`).after(`<small class='text-danger'>${msg.message}</small>`);
	setTimeout(() => {
		$(`#${msg.target}`).removeClass('error');
		$(`#${msg.target}`)
			.siblings('small')
			.remove();
	}, 5000);
});

Shiny.addCustomMessageHandler('referral-global-error', (msg) => {
	$('#referral-global-error').html(msg.message);
	setTimeout(() => {
		$('#referral-global-error').html('');
	}, 5000);
});

/* ********************* UI/BIGDASH HANDLERS **************************** */
Shiny.addCustomMessageHandler('show-tabs', (msg) => {
	setTimeout(() => {
		$('.sidebar-content')
			.children()
			.each((index, el) => {
				if($(el).hasClass('collapse'))
					return;

				$(el).show();
			});

	$('#sidebar-container .collapse')
		.first()
		.find('hr')
		.first()
		.hide();

	if(!$('.big-tab[data-name="load-tab"]').is(':visible'))
		return;

	$('.tab-trigger[data-target="dataview-tab"]').trigger('click');
	$('#sidebar-help-container').show();
	}, 1000);
});

Shiny.addCustomMessageHandler('bigdash-select-tab', (msg) => {
    $(`.tab-trigger[data-target=${msg.value}]`).trigger('click');
});

Shiny.addCustomMessageHandler('bigdash-hide-menuitem', (msg) => {
    $(`.tab-trigger[data-target=${msg.value}]`).hide();
    $(`.tab-trigger-hr[data-target=${msg.value}]`).hide();
});

Shiny.addCustomMessageHandler('bigdash-show-menuitem', (msg) => {
    $(`.tab-trigger[data-target=${msg.value}]`).show();
	$(`.tab-trigger-hr[data-target=${msg.value}]`).show();
});

Shiny.addCustomMessageHandler('bigdash-hide-tab', (msg) => {
    $(`.big-tab[data-name=${msg.value}]`).hide();
});

Shiny.addCustomMessageHandler('bigdash-show-tab', (msg) => {
    $(`.big-tab[data-name=${msg.value}]`).show();
});


/* ********************* HUBSPOT HANDLER **************************** */
$(document).ready(function() {
    /* From https://stackoverflow.com/questions/74643167 */
    $(".tab-trigger").on("click", function(e) {
	var tabId = $(e.target).data("target");
	let hasHsq = (typeof window._hsq !== 'undefined' && window._hsq !== null)
 	/* https://developers.hubspot.com/docs/api/events/tracking-code#tracking-in-single-page-applications */
//	console.log('[tab-trigger:click] user = ' + user);	
	if(hasHsq && user !== '' && user !== 'undefined') {
	    var _hsq = window._hsq = window._hsq || [];
	    var orginalTitle = document.title;
	    document.title = orginalTitle + ' > ' + tabId ;
//	    console.log('[_hsq.push] ' + document.title);	    
	    _hsq.push(["identify", { email: user }]);  // set to current user
	    _hsq.push(['setPath', '#' + tabId]);
	    _hsq.push(['trackPageView']);
	    document.title = orginalTitle;
	}
    });

});




document.addEventListener('shown.bs.popover', function(event) {
	// Get the popover container (Bootstrap 5 attaches the popover to the body by default)
	var popover = document.body.querySelector('.popover');
  
	if (!popover) {
	  console.log("Popover not found.");
	  return;
	}
  
	// Find the download button inside the popover using its class
	var targetNode = popover.querySelector('.shiny-download-link');
  
	if (!targetNode) {
	  console.log("Download button not found.");
	  return;
	}
  
	// Options for the observer (which attributes to watch)
	var config = { attributes: true, attributeFilter: ['href'] };
  
	// Callback function to execute when mutations are observed
	var callback = function(mutationsList, observer) {
	  for (var mutation of mutationsList) {
		if (mutation.type === 'attributes' && mutation.attributeName === 'href') {
		  // Enable the button by removing the "disabled" attribute
		  $(targetNode).removeClass('disabled');
		}
	  }
	};
  
	// Create an observer instance linked to the callback function
	var observer = new MutationObserver(callback);
  
	// Start observing the target node for configured mutations
	observer.observe(targetNode, config);
  });