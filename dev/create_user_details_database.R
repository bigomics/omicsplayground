# Step 1: Install and load RSQLite
library(RSQLite)

# Step 2: Create a new SQLite database and establish a connection
connection <- dbConnect(RSQLite::SQLite(), dbname = "etc/user_details.sqlite")

# Step 3: Create a table with the specified fields
dbExecute(connection, "
  CREATE TABLE IF NOT EXISTS users (
    username TEXT,
    email TEXT UNIQUE, -- email cannot be duplicated
    password TEXT,
    expiry DATE,
    level INTEGER,
    user_limit INTEGER,
    ENABLE_CHIRP INTEGER,  -- Boolean fields stored as integers
    ENABLE_DELETE INTEGER,
    ENABLE_PGX_DOWNLOAD INTEGER,
    ENABLE_PUBLIC_SHARE INTEGER,
    ENABLE_UPLOAD INTEGER,
    ENABLE_USER_SHARE INTEGER,
    MAX_DATASETS INTEGER,  -- Integer fields
    MAX_SAMPLES INTEGER,
    MAX_COMPARISONS INTEGER,
    MAX_GENES INTEGER,
    MAX_GENESETS INTEGER,
    MAX_SHARED_QUEUE INTEGER,
    TIMEOUT INTEGER,
    WATERMARK INTEGER
  )
")

# Step 4: Close the database connection
dbDisconnect(connection)

## Population example

# Re-establish the connection to the SQLite database
connection <- dbConnect(RSQLite::SQLite(), dbname = "etc/user_details.sqlite")

# Insert records into the table
dbExecute(connection, "
  INSERT INTO users (
    username, email, password, expiry, level, user_limit,
    ENABLE_CHIRP, ENABLE_DELETE, ENABLE_PGX_DOWNLOAD, ENABLE_PUBLIC_SHARE, ENABLE_UPLOAD, ENABLE_USER_SHARE,
    MAX_DATASETS, MAX_SAMPLES, MAX_COMPARISONS, MAX_GENES, MAX_GENESETS, MAX_SHARED_QUEUE, TIMEOUT, WATERMARK
  ) VALUES 
  ('user1', 'test@demo.com', 'pass123', '2022-12-31', 1, 5, 
   1, 0, 1, 0, 1, 0, 
   10, 100, 5, 5000, 50, 2, 30, 1),
  ('user2', 'user2@example.com', 'pass456', '2024-01-31', 2, 10, 
   0, 1, 0, 1, 0, 1, 
   20, 200, 10, 10000, 100, 4, 60, 0)
")

# Close the database connection
dbDisconnect(connection)