(ns tservice.plugins.quartet-protqc-report.protqc
  "A wrapper for protqc tool."
  (:require [tservice.api.config :refer [add-env-to-path]]
            [tservice.lib.files :refer [is-localpath?]]
            [clojure.string :as clj-str]
            [clojure.java.shell :as shell :refer [sh]]))

(defn call-protqc!
  "Call protqc bash script. more details on https://github.com/chinese-quartet/ProtQC
   exp-file: Proteomics profiled data. 
   meta-file: proteomics metadata.
   result-dir: A directory for result files.
  "
  [exp-file meta-file result-dir]
  (shell/with-sh-env {:PATH   (add-env-to-path "quartet-protqc-report")
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
