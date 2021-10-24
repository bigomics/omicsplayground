GOOGLE_BASE_URL <- "https://firestore.googleapis.com/v1/projects"

google_base_url <- function(){
	sprintf(
		"%s/%s/databases/(default)/documents/users",
		GOOGLE_BASE_URL,
		Sys.getenv("OMICS_GOOGLE_PROJECT")
	)
}

google_user_get <- function(key, uid){
	if(missing(key))
		stop("Missing key")
	
	if(missing(uid))
		stop("Missing uid")

	url <- sprintf(
		"%s/%s?key=%s", 
		google_base_url(), 
		URLencode(uid), 
		key
	)

	req <- httr::GET(url)
	httr::content(req)
}

google_user_create <- function(key, email, plan = "free"){
	if(missing(key))
		stop("Missing key")
	
	url <- sprintf(
		"%s?documentId=%s&key=%s", 
		google_base_url(), 
		URLencode(email),
		key
	)

	body <- list(
		fields = list(
			plan = list(
				stringValue = plan
			)
		)
	)

	req <- httr::POST(
		url, 
		body = body, 
		encode = "json",
		httr::content_type_json()
	)
	httr::content(req)	
}
