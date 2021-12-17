(ns quartet-protqc-report.core
  (:require [tservice-core.tasks.http :as http-task]
            [quartet-protqc-report.spec :as spec]
            [quartet-protqc-report.task :as task]))

(def metadata
  (http-task/make-routes "quartet-protqc-report" :ReportPlugin
                         {:method-type :post
                          :endpoint "quartet-protqc-report"
                          :summary "Generate the QC Report for Quartet Proteomics data."
                          :body-schema spec/quartet-protqc-report-params-body
                          :response-schema any?
                          :handler task/post-handler}))

(def events-init task/events-init)
