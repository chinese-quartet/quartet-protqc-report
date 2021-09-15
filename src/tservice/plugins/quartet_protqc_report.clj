(ns tservice.plugins.quartet-protqc-report
  (:require [clojure.data.json :as json]
            [clojure.spec.alpha :as s]
            [spec-tools.core :as st]
            [clojure.tools.logging :as log]
            [tservice.lib.files :as ff]
            [tservice.api.storage.fs :as fs-api]
            [tservice.lib.fs :as fs-lib]
            [tservice.vendor.multiqc :as mq]
            [tservice.plugins.quartet-protqc-report.protqc :as protqc]
            [tservice.api.task :refer [make-events-init publish-event! make-plugin-metadata create-task! update-process!]]))

;;; ------------------------------------------------ Event Specs ------------------------------------------------
(s/def ::metadata_file
  (st/spec
   {:spec                (s/and string? #(re-matches #"^[a-zA-Z0-9]+:\/\/(\/|\.\/)[a-zA-Z0-9_]+.*" %))
    :type                :string
    :description         "File path for metadata file, such as file:///xxx/xxx/metadata.csv"
    :swagger/default     nil
    :reason              "The filepath must be string."}))

(s/def ::data_file
  (st/spec
   {:spec                (s/and string? #(re-matches #"^[a-zA-Z0-9]+:\/\/(\/|\.\/)[a-zA-Z0-9_]+.*" %))
    :type                :string
    :description         "File path for proteomics profiled data, such as file:///xxx/xxx/data.csv"
    :swagger/default     nil
    :reason              "The filepath must be string."}))

(s/def ::name
  (st/spec
   {:spec                string?
    :type                :string
    :description         "The name of the report"
    :swagger/default     ""
    :reason              "Not a valid name"}))

(s/def ::description
  (st/spec
   {:spec                string?
    :type                :string
    :description         "Description of the report"
    :swagger/default     ""
    :reason              "Not a valid description."}))

(def quartet-protqc-report-params-body
  "A spec for the body parameters."
  (s/keys :req-un [::name ::data_file ::metadata_file]
          :opt-un [::description]))

;;; ------------------------------------------------ Event Metadata ------------------------------------------------
(def metadata
  (make-plugin-metadata
   {:name "quartet-protqc-report"
    :summary "Visualizes Quality Control(QC) results from proteomics data for Quartet Project."
    :params-schema quartet-protqc-report-params-body
    :handler (fn [{:keys [name data_file metadata_file description]
                   :or {description (format "Quality control report for %s" name)}
                   :as payload}]
               (let [payload (merge {:description description} payload)
                     data-file (protqc/correct-filepath data_file)
                     metadata-file (protqc/correct-filepath metadata_file)
                     workdir (fs-api/dirname data-file)
                     log-path (fs-lib/join-paths workdir "log")
                     response {:report (format "%s/multireport.html" workdir)
                               :log log-path
                               :response-type :data2report}
                     task-id (create-task! {:name           name
                                            :description    description
                                            :payload        payload
                                            :plugin-name    "quartet-protqc-report"
                                            :plugin-type    "ReportPlugin"
                                            :plugin-version "v0.1.0"
                                            :response       response})
                     result-dir (fs-lib/join-paths workdir "results")]
                 (fs-lib/create-directories! result-dir)
                 (log/info (format "Create a report %s with %s" name payload))
                 (spit log-path (json/write-str {:status "Running"
                                                 :msg ""}))
                 (update-process! task-id 0)
                 (publish-event! "quartet_protqc_report"
                                 {:data-file data-file
                                  :metadata-file metadata-file
                                  :dest-dir workdir
                                  :task-id task-id
                                  :metadata {:name name
                                             :description description
                                             :plugin-name "quartet-protqc-report"
                                             :plutin-type "ReportPlugin"
                                             :plugin-version "v0.1.0"}})
                 response))
    :plugin-type :ReportPlugin
    :response-type :data2report}))

(defn update-log-process!
  "Update message into log file and process into database."
  [log-path coll task-id process]
  (spit log-path (json/write-str coll))
  (update-process! task-id process))

(defn date
  []
  (.format (java.text.SimpleDateFormat. "yyyy-MM-dd")
           (new java.util.Date)))

(defn- make-report!
  [{:keys [data-file metadata-file dest-dir metadata task-id]}]
  (let [log-path (fs-lib/join-paths dest-dir "log")
        result-dir (fs-lib/join-paths dest-dir "results")
        parameters-file (fs-lib/join-paths result-dir "general_information.json")
        results (ff/chain-fn-coll [(fn []
                                     (update-process! task-id 20)
                                     (protqc/call-protqc! data-file metadata-file result-dir))
                                   (fn []
                                     (update-process! task-id 50)
                                     (spit parameters-file (json/write-str {"Report Name" (:name metadata)
                                                                            "Description" (:description metadata)
                                                                            "Report Tool" (format "%s-%s"
                                                                                                  (:plugin-name metadata)
                                                                                                  (:plugin-version metadata))
                                                                            "Team" "Quartet Team"
                                                                            "Date" (date)}))
                                     {:status "Success" :msg ""})
                                   (fn []
                                     (update-process! task-id 80)
                                     (mq/multiqc result-dir dest-dir {:template "default" :title "Quartet Report for Proteomics"}))]
                                  (fn [result] (= (:status result) "Success")))
        status (:status (last results))
        msg (apply str (map :msg results))
        process (if (= status "Success") 100 -1)]
    (log/info (format "Running batch command: %s" (pr-str results)))
    (update-log-process! log-path {:status status
                                   :msg msg}
                         task-id process)))

;;; --------------------------------------------------- Lifecycle ----------------------------------------------------

(def events-init
  "Automatically called during startup; start event listener for quartet_protqc_report events."
  (make-events-init "quartet_protqc_report" make-report!))
