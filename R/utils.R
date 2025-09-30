#wrapper for the sample function that will take a sample of size 1 if x has length 1
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Clear notes in check() log
id <- t1 <- t2 <- t.first <- event.first <- NULL
