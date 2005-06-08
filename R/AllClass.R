## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("lmList.confint", contains = "array")



