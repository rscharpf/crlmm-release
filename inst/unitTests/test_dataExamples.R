test_dataExamples <- function(){
	data(cnSetExample)
	checkTrue(validObject(cnSetExample))
	data(cnSetExample2)
	checkTrue(validObject(cnSetExample2))
}
