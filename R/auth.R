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

google_user_get <- function(credentials, email){
	if(missing(credentials))
		stop("Missing credentials")
	
	url <- sprintf("%s/%s", google_base_url(), URLencode(email))

	req <- httr::GET(
		url,
		httr::config(token = credentials)
	)
	httr::content(req)
}

google_user_patch <- function(credentials, user){
	email <- URLencode(user$email)
	url <- sprintf("%s&document.name=%s", google_base_url(), email)
}
