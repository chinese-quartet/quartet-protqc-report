(ns quartet-protqc-report.protqc
  "A wrapper for protqc tool."
  (:require [clojure.string :as clj-str]
            [local-fs.core :as fs-lib]
            [clojure.java.shell :as shell :refer [sh]]
            [tservice-core.plugins.util :refer [call-command!]]
            [tservice-core.plugins.env :refer [get-context-path add-env-to-path]]
            [clojure.tools.logging :as log]
            [quartet-protqc-report.version :as v]))

(defn call-protqc!
  "Call protqc bash script. more details on https://github.com/chinese-quartet/ProtQC
   exp-file: Proteomics profiled data. 
   meta-file: proteomics metadata.
   result-dir: A directory for result files."
  [exp-file meta-file result-dir]
  (let [command ["bash" "-c"
                 (format "protqc.sh -d %s -m %s -o %s" exp-file meta-file result-dir)]
        rprofile (fs-lib/join-paths (get-context-path :env v/plugin-name) "Rprofile")
        path (add-env-to-path v/plugin-name)
        ;; When you are in local mode, the context-path doesn't exist.
        rprofile (if (= rprofile (format "%s/Rprofile" v/plugin-name))
                   (System/getProperty "R_PROFILE_USER")
                   rprofile)
        path (if (= path (format "%s/bin" v/plugin-name))
               (System/getenv "PATH")
               path)]
    (log/info "PATH variable: " path)
    (log/info "Rprofile file is in " rprofile)
    (shell/with-sh-env {:PATH   path
                        :R_PROFILE_USER rprofile
                        :LC_ALL "en_US.utf-8"
                        :LANG   "en_US.utf-8"}
      (let [result (apply sh command)]
        {:status (if (= (:exit result) 0) "Success" "Error")
         :msg (str (:out result) "\n" (:err result))}))))

(defn multiqc
  "A multiqc wrapper for generating multiqc report:
   TODO: set the absolute path of multiqc binary instead of environment variable

  Required:
  analysis-dir: Analysis directory, e.g. data directory from project
  outdir: Create report in the specified output directory.

  Options:
  | key                | description |
  | -------------------|-------------|
  | :dry-run?          | Dry run mode |
  | :filename          | Report filename. Use 'stdout' to print to standard out. |
  | :comment           | Custom comment, will be printed at the top of the report. |
  | :title             | Report title. Printed as page header, used for filename if not otherwise specified. |
  | :force?            | Overwrite any existing reports |
  | :prepend-dirs?     | Prepend directory to sample names |
  | :template          | default, other custom template    |
  | :config            | Where is the config file          |
  | :env               | An environemnt map for running multiqc, such as {:PATH (get-path-variable)} |

  Example:
  (multiqc 'XXX' 'YYY' {:filename       'ZZZ'
                        :comment        ''
                        :title          ''
                        :force?         true
                        :prepend-dirs?  true})"
  [analysis-dir outdir {:keys [dry-run? filename comment title force? prepend-dirs? template config env]
                        :or   {dry-run?      false
                               force?        true
                               prepend-dirs? false
                               filename      "multiqc_report.html"
                               comment       ""
                               template      "default"
                               title         "iSEQ Analyzer Report"}}]
  (let [force-arg   (if force? "--force" "")
        dirs-arg    (if prepend-dirs? "--dirs" "")
        config-arg  (if config (str "-c " config) "")
        multiqc-command (filter #(> (count %) 0) ["multiqc"
                                                  force-arg dirs-arg config-arg
                                                  "--title" (format "'%s'" title)
                                                  "--comment" (format "'%s'" comment)
                                                  "--filename" filename
                                                  "--outdir" outdir
                                                  "-t" template
                                                  analysis-dir])
        command (clj-str/join " " multiqc-command)]
    (if dry-run?
      (log/info command)
      (if env
        (call-command! command env)
        (call-command! command)))))

(defn is-localpath?
  [filepath]
  (re-matches #"^file:\/\/.*" filepath))

(defn correct-filepath
  [filepath]
  (if (is-localpath? filepath)
    (clj-str/replace filepath #"^file:\/\/" "")
    filepath))
