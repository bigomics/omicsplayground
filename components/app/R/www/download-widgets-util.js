Shiny.addCustomMessageHandler("export-widget", function(message) {
    var widgetId = message.id;
    var plotName = message.name;
    var widgetElement = document.getElementById(widgetId);
    var scaleWidth = message.width/2;
    var scaleHeight = message.height/2;
  
    if (widgetElement) {
      domtoimage.toPng(widgetElement, {
        width: widgetElement.clientWidth * scaleWidth,
        height: widgetElement.clientHeight * scaleHeight,
        style: {
          transform: "scale(" + scaleWidth + ", " + scaleHeight + ")",
          transformOrigin: "top left"
        }
      })
      .then(function(dataUrl) {
        var link = document.createElement('a');
        link.href = dataUrl;
        link.download = plotName;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      })
      .catch(function(error) {
        console.error("Failed to export widget:", error);
      });
    } else {
      console.error("Element with ID " + widgetId + " not found.");
    }
});
