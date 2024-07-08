function copyPlotModuleInfo() {
    var elements = document.getElementsByClassName('popover-body')[0].children[0].children;
    var textToCopy = '';
    for (var i = 0; i < elements.length; i++) {
        if (elements[i].className === "plotmodule-info") {
            textToCopy += elements[i].innerText + '\n';
        }
    }
    var tempInput = document.createElement('input');
    tempInput.value = textToCopy.trim();
    document.body.appendChild(tempInput);
    tempInput.select();
    document.execCommand('copy');
    document.body.removeChild(tempInput);
}
