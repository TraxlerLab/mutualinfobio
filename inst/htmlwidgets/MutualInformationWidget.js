HTMLWidgets.widget({

  name: 'MutualInformationWidget',

  type: 'output',

  factory: function(el, width, height) {


    return {

      renderValue: function(x) {

        printJSAV('el', x.data, x.options);

      },

      resize: function(width, height) {

        // TODO: code to re-render the widget with a new size

      }

    };
  }
});
