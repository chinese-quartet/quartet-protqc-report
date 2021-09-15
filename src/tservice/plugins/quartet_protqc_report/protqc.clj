(ns tservice.plugins.quartet-protqc-report.protqc
  "A wrapper for protqc tool."
  (:require [tservice.api.config :refer [add-env-to-path]]
            [tservice.lib.files :refer [is-localpath? get-plugin-jar-env-dir]]
            [clojure.string :as clj-str]
            [tservice.lib.fs :as fs-lib]
            [clojure.java.shell :as shell :refer [sh]]
            [clojure.java.io :refer [file]]))

(defn call-protqc!
  "Call protqc bash script. more details on https://github.com/chinese-quartet/ProtQC
   exp-file: Proteomics profiled data. 
   meta-file: proteomics metadata.
   result-dir: A directory for result files.
  "
  [exp-file meta-file result-dir]
  (shell/with-sh-env {:PATH   (add-env-to-path "quartet-protqc-report")
                      :R_PROFILE_USER (fs-lib/join-paths (get-plugin-jar-env-dir "quartet-protqc-report") "Rprofile")
                      :LC_ALL "en_US.utf-8"
                      :LANG   "en_US.utf-8"}
    (let [command ["bash" "-c"
                   (format "protqc.sh -d %s -m %s -o %s" exp-file meta-file result-dir)]
          result  (apply sh command)
          status (if (= (:exit result) 0) "Success" "Error")
          msg (str (:out result) "\n" (:err result))]
      {:status status
       :msg msg})))

(defn correct-filepath
  [filepath]
  (if (is-localpath? filepath)
    (clj-str/replace filepath #"^file:\/\/" "")
    filepath))

(defn ^String basename
  "Returns the basename of 'path'.
   This works by calling getName() on a java.io.File instance. It's prefered
   over last-dir-in-path for that reason.
   Parameters:
     path - String containing the path for an item in iRODS.
   Returns:
     String containing the basename of path."
  [^String path]
  (.getName (file path)))

(defn ^String dirname
  "Returns the dirname of 'path'.
   This works by calling getParent() on a java.io.File instance.
   Parameters:
     path - String containing the path for an item in iRODS.
   Returns:
     String containing the dirname of path."
  [^String path]
  (when path (.getParent (file path))))
