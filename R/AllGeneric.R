#if (!isGeneric("getGroups")) {
#    ## Return the groups associated with object according to form.
#    setGeneric("getGroups",
#               function(object, form, level, data, sep, ...)
#               standardGeneric("getGroups"))
#}

#if (!isGeneric("getGroupsFormula")) {
#    ## Return the formula(s) for the groups associated with object.
#    ## The result is a one-sided formula unless asList is TRUE in which case
#    ## it is a list of formulas, one for each level.
#    setGeneric("getGroupsFormula",
#               function(object, asList = FALSE, sep = "/")
#               standardGeneric("getGroupsFormula"))
#}

#if (!isGeneric("getCovariate")) {
#    ## Return the primary covariate associated with object
#    setGeneric("getCovariate",
#               function(object, form = formula(object), data = list())
#               standardGeneric("getCovariate"))
#}

#if (!isGeneric("getResponse")) {
#    ## Return the response variable associated with object
#    setGeneric("getResponse",
#               function(object, form = formula(object))
#               standardGeneric("getResponse"))
#}


if (!isGeneric("lmList")) {
    setGeneric("lmList",
               function(formula, data, family, subset, weights,
                        na.action, offset, pool, ...)
               standardGeneric("lmList"))
}

if (!isGeneric("pooledSD")) {
    setGeneric("pooledSD", function(object) standardGeneric("pooledSD"))
}

if (!isGeneric("gsummary")) {
    setGeneric("gsummary",
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
               standardGeneric("gsummary"))
}

