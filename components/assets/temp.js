let db;
let pricing;
let user;
var dimension = [0, 0];

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

	// Two non-defaults, both so a data-bs-placement="left" tooltip in the 12rem
	// settings panel actually opens on the left. container: 'body' keeps the
	// panel from being the clipping boundary, and flip's altAxis check goes off
	// because these tooltips are long paragraphs: at the 200px the stylesheet
	// allows them they run several hundred pixels tall, so one centred on an
	// input near the top of the panel pokes out above the viewport. flip counts
	// that overflow - perpendicular to the side asked for, with 1300px of room
	// still free on the left - as "left does not fit", walks the fallbacks and
	// lands on bottom. With the check off the placement survives and
	// preventOverflow just slides the tooltip down into view.
	const popperConfig = (defaults) => ({
		...defaults,
		modifiers: [...defaults.modifiers, { name: 'flip', options: { altAxis: false } }]
	});
	const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
	const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl, { container: 'body', popperConfig }));
})


Shiny.addCustomMessageHandler('enableInfo', (data) => {
  Shiny.setInputValue(data.id, data.value);
});

Shiny.addCustomMessageHandler('copilot-chat-batch-append', (msg) => {
  const chatId = msg && msg.id ? msg.id : null;
  const chatEl = chatId ? document.getElementById(chatId) : null;
  const messages = Array.isArray(msg && msg.messages) ? msg.messages : [];
  const doneInputId = msg && msg.done_input_id ? msg.done_input_id : null;
  const rawBatchSize = msg && Number.isFinite(msg.batch_size) ? Number(msg.batch_size) : 32;
  const batchSize = Math.max(1, rawBatchSize);

  if (!chatEl || messages.length === 0) {
    if (doneInputId) {
      Shiny.setInputValue(doneInputId, Date.now(), { priority: 'event' });
    }
    return;
  }

  let idx = 0;
  const appendNext = () => {
    const limit = Math.min(idx + batchSize, messages.length);
    for (; idx < limit; idx++) {
      const message = messages[idx] || {};
      const detail = {
        role: message.role === 'user' ? 'user' : 'assistant',
        content: typeof message.content === 'string' ? message.content : '',
        content_type: 'markdown',
        html_deps: '[]',
        chunk_type: null,
        operation: null
      };
      const ev = new CustomEvent('shiny-chat-append-message', { detail });
      chatEl.dispatchEvent(ev);
    }

    if (idx < messages.length) {
      window.requestAnimationFrame(appendNext);
      return;
    }

    if (doneInputId) {
      Shiny.setInputValue(doneInputId, Date.now(), { priority: 'event' });
    }
  };

  appendNext();
});

const logout = () => {
	unloadSidebar();
	sidebarClose();
 	Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

const logoutInApp = () => {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", "cookie_remove", true);
	xhr.onload = function() {
		window.location = window.location.href;
    };
	xhr.send();
};

const quit = () => {
        Shiny.setInputValue('quit', 1, {priority: 'event'});
};

Shiny.addCustomMessageHandler('shinyproxy-logout', (msg) => {
        window.location.assign("/logout");
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

Shiny.addCustomMessageHandler('redirect', function(message) {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", "cookie", true);
	xhr.setRequestHeader("Header-User-Cookie", message);
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200)
		window.location = message;
	};
xhr.send();
});

Shiny.addCustomMessageHandler('redirect_nonce', function(message) {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", "cookie_nonce", true);
	xhr.setRequestHeader("Header-User-Cookie", message);
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200)
		window.location = message;
	};
	xhr.send();
});

var hrefUpdatedStatus = {};  // Object to track href update status for each popover

document.addEventListener('shown.bs.popover', function(event) {
	// The popover that was just shown, found through the trigger's own
	// aria-describedby (Bootstrap points it at the tip while that tip is open).
	// document.body.querySelector('.popover') walked the whole page instead --
	// every board and every node inside every drawn SVG -- to reach a popover
	// Bootstrap appends at the end of <body>, and returned the wrong one
	// whenever more than one was open.
	var trigger = event.target;
	var popoverId = trigger && trigger.getAttribute && trigger.getAttribute('aria-describedby');
	var popover = popoverId ? document.getElementById(popoverId) : null;

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

	// Get the id attribute of the download button
    var buttonId = targetNode.getAttribute('id');
    if (!buttonId) {
        console.log("Download button ID not found.");
        return;
    }

	// If href has been updated previously for this popover, always remove the disabled class
    if (hrefUpdatedStatus[buttonId]) {
        //$(targetNode).removeClass('disabled');
		return;
    } else {
		$(targetNode).addClass('disabled');
	}
  
	// Options for the observer (which attributes to watch)
	var config = { attributes: true, attributeFilter: ['href'] };
  
	// Callback function to execute when mutations are observed
	var callback = function(mutationsList, observer) {
	  for (var mutation of mutationsList) {
		if (mutation.type === 'attributes' && mutation.attributeName === 'href') {
		  // Enable the button by removing the "disabled" attribute
		  $(targetNode).removeClass('disabled');
		  // Set the flag to true as href has been updated for this popover
		  hrefUpdatedStatus[buttonId] = true;
		}
	  }
	};
  
	// Create an observer instance linked to the callback function
	var observer = new MutationObserver(callback);
  
	// Start observing the target node for configured mutations
	observer.observe(targetNode, config);

	// Dropped when the popover closes. Without this every single show left a
	// live observer behind, for the whole session.
	trigger.addEventListener('hidden.bs.popover', function () {
	  observer.disconnect();
	}, { once: true });
  });
