Shiny.addCustomMessageHandler("export-widget", function(message) {
    var widgetId = message.id;
    var plotName = message.name;
    var widgetElement = document.getElementById(widgetId);
    var scaleWidth = message.width/2;
    var scaleHeight = message.height/2;
    var format = message.format || "png";

    if (widgetElement) {
      if (format === "png") {
          // Export as PNG
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
              link.download = plotName + ".png";
              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
          })
          .catch(function(error) {
              console.error("Failed to export widget as PNG:", error);
          });
      } else if (format === "pdf") {
          // Export as PDF
          domtoimage.toPng(widgetElement, {
              width: widgetElement.clientWidth * scaleWidth,
              height: widgetElement.clientHeight * scaleHeight,
              style: {
                  transform: "scale(" + scaleWidth + ", " + scaleHeight + ")",
                  transformOrigin: "top left"
              }
          })
          .then(function(dataUrl) {
              const { jsPDF } = window.jspdf;
              var pdf = new jsPDF({
                  orientation: 'landscape',
                  unit: 'px',
                  format: [widgetElement.clientWidth, widgetElement.clientHeight]
              });

              pdf.addImage(dataUrl, 'PNG', 0, 0, widgetElement.clientWidth, widgetElement.clientHeight);
              pdf.save(plotName + ".pdf");
          })
          .catch(function(error) {
              console.error("Failed to export widget as PDF:", error);
          });
      }
  } else {
      console.error("Element with ID " + widgetId + " not found.");
  }
});
