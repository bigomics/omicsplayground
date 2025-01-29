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


const unloadSidebar = () => {
	$('.sidebar-content')
		.children()
		.each((index, el) => {
			if(index == 0 || index == 1){
				$(el).show();
				return;
			}
			if($(el).hasClass('collapse')){
				$(el).removeClass('show');
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

const logout = () => {
	unloadSidebar();
	sidebarClose();
 	Shiny.setInputValue('userLogout', 1, {priority: 'event'});
};

const logoutInApp = () => {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", "/cookie_remove", true);
	xhr.onload = function() {
		window.location = window.location.origin + '/';
    };
	xhr.send();
};

const quit = () => {
        Shiny.setInputValue('quit', 1, {priority: 'event'});
};

Shiny.addCustomMessageHandler('shinyproxy-logout', (msg) => {
        window.location.assign("/logout");
});

/* ********************* UI/BIGDASH HANDLERS **************************** */
Shiny.addCustomMessageHandler('show-tabs', (msg) => {
	setTimeout(() => {
	$('.sidebar-content')
		.children()
		.each((index, el) => {
      if ($(el).hasClass('collapse')) { // && !$(el).hasClass('nodisp')) {
			// if($(el).hasClass('collapse')) {
				// if($(el).hasClass('show'))
				// 	$(el).removeClass('show');
				//$(el).hide();
				$(el).removeClass('show');
				$(el).css({'display' : ''});
				return;
			}
			if($(el).hasClass('w-100')) {
				// console.log($(el).children().children()[[1]]);
				$(el).children().children()[[1]].classList.remove('fa-angle-down');
				$(el).children().children()[[1]].classList.add('fa-angle-right');
			}
      if (!$(el).hasClass('nodisp')) {
        $(el).show();
      }
		});

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

Shiny.addCustomMessageHandler('bigdash-remove-tab', (msg) => {
    $(`.big-tab[data-name=${msg.value}]`).remove();
	$(`.tab-settings:has(a#${msg.value.slice(0, -4)}-options)`).remove();
    // $(`[data-target=${msg.value}]`).addClass("nodisp");
    $(`[data-target=${msg.value}]`).hide();
});

Shiny.addCustomMessageHandler('bigdash-hide-menu-element', (msg) => {
    $(`span:contains(${msg.value})`).closest('p').hide();
    $(`span:contains(${msg.value})`).closest('p').addClass("nodisp");
});

Shiny.addCustomMessageHandler('bigdash-show-menu-element', (msg) => {
    $(`span:contains(${msg.value})`).closest('p').show();
    $(`span:contains(${msg.value})`).closest('p').removeClass("nodisp");
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
  });
