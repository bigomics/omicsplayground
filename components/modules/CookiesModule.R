# Get desired cookie from all the available ones
extract_cookie_value <- function(session, cookie_name) {
  cookie_string <- session$request$HTTP_COOKIE
  if (is.null(cookie_string)) {
    return(NULL)
  }
  cookies <- unlist(strsplit(cookie_string, "; "))
  user_cookie <- grep(paste0("^", cookie_name, "="), cookies, value = TRUE)
  if (length(user_cookie) > 0) {
    return(unlist(strsplit(user_cookie, "="))[2])
  } else {
    return(NULL)
  }
}

# Decrypt cookie
decrypt_cookie <- function(cookie, nonce) {
  key_base64 <- readLines(file.path(OPG, "etc/keys/cookie.txt"))[1]
  email_nonce_raw <- tryCatch(
    {
      sodium::hex2bin(nonce)
    },
    error = function(w) {
      ""
    }
  )
  email_raw <- tryCatch(
    {
      sodium::hex2bin(cookie)
    },
    error = function(w) {
      ""
    }
  )
  attr(email_raw, "nonce") <- email_nonce_raw
  # Return NULL if not successful decryption
  email <- tryCatch(
    {
      unserialize(
        sodium::data_decrypt(
          email_raw,
          key = sodium::sha256(charToRaw(key_base64))
        )
      )
    },
    error = function(w) {
      NULL
    }
  )
}

# Get cookies and decrypt them
get_and_decrypt_cookie <- function(session) {
  cookie <- extract_cookie_value(session, "persistentOPG")
  nonce <- extract_cookie_value(session, "persistentOPG_nonce")
  if (!is.null(cookie) & !is.null(nonce)) {
    decrypted_cookie <- decrypt_cookie(cookie, nonce)
    return(decrypted_cookie)
  } else {
    return(NULL)
  }
}

# Save encrypted session cookie
save_session_cookie <- function(session, cred) {
  key_base64 <- readLines(paste0(OPG, "/etc/keys/cookie.txt"))[1]
  passkey <- sodium::sha256(charToRaw(key_base64))
  plaintext <- isolate(cred$email)
  plaintext.raw <- serialize(plaintext, NULL)
  ciphertext <- sodium::data_encrypt(plaintext.raw, key = passkey)

  email_encrypted_nonce <- paste(as.character(attr(ciphertext, "nonce")), collapse = "")
  email_encrypted_value <- paste(as.character(ciphertext), collapse = "")

  session$sendCustomMessage(type = "redirect", message = email_encrypted_value)
  session$sendCustomMessage(type = "redirect_nonce", message = email_encrypted_nonce)
}
