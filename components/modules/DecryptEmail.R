decrypt_email <- function(mail, nonce) {
  key_base64 <- readLines(file.path(OPG, "etc/keys/email.txt"))[1]
  if(is.null(mail) | is.null(nonce)) {
    return(NULL)
  }
  email_nonce_raw <- sodium::hex2bin(nonce)
  email_raw <- sodium::hex2bin(mail)
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
  return(email)
}