# Block R's download, URL and installer entry points during a cached mapper trial.
for (name in c("download.file", "install.packages")) {
  trace(name, where = asNamespace("utils"),
        tracer = quote(stop("Offline audit: network/install entry point blocked", call. = FALSE)),
        print = FALSE)
}
trace("url", where = baseenv(),
      tracer = quote(stop("Offline audit: URL entry point blocked", call. = FALSE)),
      print = FALSE)
source(Sys.getenv("FAPROTAX2_MAPPER"), local = FALSE)
