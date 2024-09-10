Shiny.addCustomMessageHandler("export-widget", function(message) {
    var widgetId = message.id;
    var plotName = message.name;
    var widgetElement = document.getElementById(widgetId);
  
    if (widgetElement) {
      html2canvas(widgetElement).then(canvas => {
        var imgData = canvas.toDataURL("image/png");
        var link = document.createElement('a');
        link.href = imgData;
        link.download = plotName;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    });
    } else {
      console.error("Element with ID " + widgetId + " not found.");
    }
  })
