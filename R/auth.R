GOOGLE_BASE_URL <- "https://firestore.googleapis.com/v1/projects"

google_base_url <- function(){
	sprintf(
		"%s/%s/databases/(default)/documents/users",
		GOOGLE_BASE_URL,
		Sys.getenv("OMICS_GOOGLE_PROJECT")
	)
}

google_oauth <- function(){
	myapp <- httr::oauth_app("google",
		key = Sys.getenv("OMICS_GOOGLE_KEY"),
		secret = Sys.getenv("OMICS_GOOGLE_SECRET")
	)

	httr::oauth2.0_token(
		httr::oauth_endpoints("google"),
		cache = TRUE, 
		myapp,
		scope = "https://www.googleapis.com/auth/datastore"
	)
}

google_user_get <- function(credentials, user){
	if(missing(credentials))
		stop("Missing credentials")
	
	url <- sprintf("%s/%s", google_base_url(), user)

	req <- httr::GET(
		url,
		httr::config(token = credentials)
	)
	httr::content(req)
}
