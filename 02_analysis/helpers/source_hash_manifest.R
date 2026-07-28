## Verify and record hashes for cross-compartment source artifacts.

source_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("Package 'digest' is required for source hash checks.", call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

manifest_key <- function(path, root = getwd()) {
  np <- normalizePath(path, mustWork = FALSE)
  nr <- normalizePath(root, mustWork = FALSE)
  prefix <- paste0(nr, .Platform$file.sep)
  if (startsWith(np, prefix)) sub(prefix, "", np, fixed = TRUE) else path
}

verify_source_hash <- function(path, source_label, manifest_path, root = getwd()) {
  if (!file.exists(path))
    stop("cross-compartment source absent: ", path, call. = FALSE)
  if (!file.exists(manifest_path))
    stop("source hash manifest absent: ", manifest_path,
         ". Pin the cross-compartment source before reading it.", call. = FALSE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  source_path <- manifest_key(path, root = root)
  hit <- manifest[manifest$source_label == source_label &
                    manifest$source_path == source_path, , drop = FALSE]
  if (nrow(hit) != 1L)
    stop("source hash pin missing or duplicated for ", source_label, " (",
         source_path, ") in ", manifest_path, call. = FALSE)
  current <- source_sha256(path)
  expected <- hit$sha256[1]
  if (!identical(current, expected))
    stop("source hash mismatch for ", source_label, ": ", source_path,
         "\nexpected ", expected, "\nobserved ", current,
         "\nRegenerate or update the consuming stage only after reviewing the source change.",
         call. = FALSE)
  current
}

verify_optional_source_hash <- function(path, source_label, manifest_path, root = getwd()) {
  source_path <- manifest_key(path, root = root)
  if (!file.exists(path)) {
    return(data.frame(
      source_label = source_label,
      source_path = source_path,
      status = "missing",
      sha256 = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  hash <- verify_source_hash(path, source_label, manifest_path, root = root)
  data.frame(
    source_label = source_label,
    source_path = source_path,
    status = "read",
    sha256 = hash,
    stringsAsFactors = FALSE
  )
}
