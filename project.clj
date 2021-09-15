(defproject tservice/quartet-protqc-report "v0.1.0"
  :description "Visualizes Quality Control(QC) results for Quartet Project."
  :url "https://github.com/tservice-plugins/quartet-protqc-report"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :min-lein-version "2.5.0"
  :deployable false

  :dependencies
  [[org.clojure/data.csv "1.0.0"]
   [me.raynes/fs "1.4.6"]
   [org.clojure/tools.logging "1.1.0"]
   [org.clojure/core.async "0.4.500"
    :exclusions [org.clojure/tools.reader]]]

  :profiles
  {:provided
   {:dependencies
    [[org.clojure/clojure "1.10.1"]
     [org.clojars.yjcyxky/tservice "0.5.8"]]}

   :uberjar
   {:auto-clean    true
    :aot           :all
    :omit-source   true
    :javac-options ["-target" "1.8", "-source" "1.8"]
    :target-path   "target/%s"
    :resource-paths ["resources"]
    :uberjar-name  "quartet-protqc-report.tservice-plugin.jar"}})
