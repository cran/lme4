# Provisional method.
# The definition of the method should change when groupedData is defined
#  as an S4 class.

setMethod("formula",
          "groupedData",
          function(x, ...) attr(x, "formula"),
          valueClass = "formula")
