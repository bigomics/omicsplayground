// This snippet ensures that when clicking outside a PlotModule window
// if there is any dropdown active, the `active` class will be removed,
// therefore it will look greyed out again

document.body.addEventListener('click', function (evt) {
    let arr2 = ['fa-info', 'fa-bars', 'fa-download', 'fa-window-maximize', 'dropdown-button'];
    let arr1 = Array.from(evt.srcElement.classList)
    if(!arr1.some(r=> arr2.includes(r))){
       $('.show.dropdown-button.active').toggleClass('active');
    }
});

// This JS code prevents that click *inside* the dropdown itself to close it.
// It improves user experience as without this, trying to click and select the info text
// would close the dropdown. Also when using a shiny::selectInput inside the dropdown, it
// was unusable as trying to interact with it would close the dropdown.
$('.dropdown-menu').click(function(e) {
  e.stopPropagation();
});