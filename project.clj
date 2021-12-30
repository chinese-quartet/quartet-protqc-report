(defproject quartet-protqc-report "0.1.3"
  :description "Visualizes Quality Control(QC) results for Quartet Project."
  :url "https://github.com/chinese-quartet/quartet-protqc-report"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :min-lein-version "2.5.0"
  :deployable false

  :dependencies
  [[org.clojure/data.csv "1.0.0"]
   [com.github.yjcyxky/local-fs "0.1.5"]
   [org.clojure/tools.logging "1.1.0"]
   [org.clojure/tools.cli "1.0.194"]
   [metosin/spec-tools "0.10.5"]
   [org.clojure/spec.alpha "0.3.214"]]

  :plugins [[lein-cloverage "1.0.13"]
            [lein-shell "0.5.0"]
            [lein-changelog "0.3.2"]]

  :aliases {"update-version" ["shell" "sed" "-i" "" "s/version \"[0-9.]*\"/version \"${:version}\"/" "src/quartet_protqc_report/version.clj"]}
  :deploy-repositories [["releases" :clojars]]
  :release-tasks [["change" "version" "leiningen.release/bump-version"]
                  ["change" "version" "leiningen.release/bump-version" "release"]
                  ["changelog" "release"]
                  ["update-version"]
                  ["deploy"]]

  :main ^:skip-aot quartet-protqc-report.cli
  :target-path   "target/%s"
  :resource-paths ["resources"]
  :source-paths ["src"]
  :test-paths ["test"]

  :profiles
  {:provided
   {:dependencies
    [[org.clojure/clojure "1.10.1"]
     [com.github.yjcyxky/tservice-core "0.2.0"]]}

   :uberjar
   {:dependencies
    [[org.clojure/clojure "1.10.1"]
     [com.github.yjcyxky/tservice-core "0.2.0"]]
    :aot           :all
    :omit-source   false
    :javac-options ["-target" "1.8", "-source" "1.8"]
    :target-path   "target/%s"
    :resource-paths ["resources"]}})
