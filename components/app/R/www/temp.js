let db;
let pricing;
let user;

$(document).on('shiny:connected', function() {

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

    Shiny.addCustomMessageHandler('window_resize', function(message) {
        if (message.resize) $(window).resize();
    })

    Shiny.addCustomMessageHandler('set-user', (msg) => {
        $('#authentication-user').text(msg.user);
        user = msg.user;
	pricing = msg.pricing;
    });

});  // end of on.shiny.connected

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
	}, 250);

	const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
	const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));
})

Shiny.addCustomMessageHandler('manage-sub', (msg) => {
	window.location.assign(msg);
});

const logout = () => {
	Shiny.setInputValue('auth-userLogout', 1, {priority: 'event'});
	Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

const logoutInApp = () => {
    $(".tab-trigger[data-target='welcome-tab']").trigger('click');
	Shiny.setInputValue('auth-userLogout', 1, {priority: 'event'});
	Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

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

const hideSub = () => {
    $('#logSub').hide();
    $('#logSubbed').show();
    $('#sever-reload-btn').show();
}

const sendLog = () => {
  let msg  = $('#logMsg').val();

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

const sendLog2 = (msg) => {
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

$(document).ready(function() {
    console.log('*** setting up HSQ ***');
    /* Default installation */
    /* From https://stackoverflow.com/questions/74643167 */
    /*    $("a[data-toggle='tab']").on("shown.bs.tab", function(e) {*/
    $(".tab-trigger").on("click", function(e) {
	var tabId = $(e.target).data("target");
	let hasHsq = (typeof window._hsq !== 'undefined' && window._hsq !== null)
 	/* https://developers.hubspot.com/docs/api/events/tracking-code#tracking-in-single-page-applications */
	console.log('[tab-trigger:click] tabId =' + tabId);
	if(hasHsq && user !== '' && user !== 'undefined') {
	    var _hsq = window._hsq = window._hsq || [];
	    var orginalTitle = document.title;
	    document.title = orginalTitle + ' > ' + tabId ;
	    _hsq.push(["identify", { email: user }]);  // set to current user
	    _hsq.push(['setPath', '#' + tabId]);
	    _hsq.push(['trackPageView']);
	    document.title = orginalTitle;
	}
    });
});