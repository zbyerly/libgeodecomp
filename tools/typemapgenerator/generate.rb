#!/usr/bin/ruby
require 'fileutils'
require 'pathname'
require 'optparse'
require 'pp'
require 'shell'

basedir = Pathname.new($0).dirname
$: << basedir
require 'typemapgenerator'

options = {}
options[:exclude] = []
options[:extension] = "cpp"

opts = OptionParser.new do |o|
  o.banner = "Usage: #$0 [OPTIONS] PATH_TO_XML_DOC [OUTPUT_PATH]"
  o.separator "Tool for automatically generating MPI typemaps from Doxygen XML output."
  o.separator ""
  o.on("-c", "--cache CACHE_FILE",
       "Scan headers for serialization candidates and ",
       "list them in a cache file, useful for integration ",
       "into larger builds.") do |cache_file|
    options[:cache] = cache_file
  end
  o.on("-s", "--scan-only",
       "Only write out the header files containing ",
       "classes for which MPI typemaps should be ",
       "generated. May be useful for build system ",
       "integration.") do
    options[:scan] = true
  end
  o.on("-S", "--sloppy",
       "Do not perform strict member mapping, but",
       "exclude all member variables from the MPI",
       "datatype, for which either no MPI datatype",
       "is available or can be created automatically.") do
    options[:sloppy] = true
  end
  o.on("-h", "--help", "This help message") do
    puts o
    exit
  end
  o.on("-e", "--extension SUFFIX",
       "Use SUFFIX as a extension for the Typemaps",
       "source file. Defaults to \"cpp\", leading",
       "to \"typemaps.cpp\" being created.") do |suffix|
    options[:extension] = suffix
  end
  o.on("-n", "--namespace <ns>",
       "Encapsulate the Typemaps class in namespace <ns>",
       "and look after classes in <ns> whilst resolving references to them.") do |namespace|
    options[:namespace] = namespace
  end
  o.on("--header-fix REPLACEMENT_PATTERN",
       "replace a pattern in the header paths") do |pattern|
    pattern =~ /(.+):(.+)/
    options[:header_pattern] = $1
    options[:header_replacement] = $2
  end
  o.on("--macro-guard MACRO",
       "encapsulate code in #ifdef(MACRO), #endif guards") do |macro|
    options[:macro_guard] = macro
  end
  o.on("-p", "--profile",
       "profile execution") do
    options[:profiling] = true
  end
end

opts.parse!(ARGV)

if ARGV.size < 1
  STDERR.puts "#$0: PATH_TO_XML_DOC expected"
  STDERR.puts opts
  exit 1
end

xml_path = ARGV[0]
output_path = Pathname.new(ARGV[1] || "./")

if options[:cache]
  cache_file = options[:cache]
  source_dir = ARGV[0]
  headers = ARGV[1].split(":")
  FileUtils::cd(source_dir)

  typemaps_deps = []
  if File.exists?(cache_file)
    typemaps_deps = File.read(cache_file).split("\n")
  end

  # puts "headers:"
  # pp headers
  # puts "typemaps_deps:"
  # pp typemaps_deps

  # at least create a empty file to indicate that there are no
  # classes to be serialized yet...
  unless File.exists?(cache_file)
    FileUtils.touch cache_file
    # ...but fake the modification time so all headers are considered newer.
    File.utime(Time.at(0), Time.at(0), cache_file)
  end

  cache_age = File.mtime(cache_file)
  newer_headers = headers.find_all { |h| File.mtime(h) > cache_age }

  buffer = []
  if newer_headers.size > 0
    buffer = `grep "^.*friend class Typemaps;" #{newer_headers.join(" ")} 2>/dev/null`

    buffer = buffer.split("\n").map do |h|
      h =~ /^(.+)\: /
      $1
    end
  end

  if buffer.size > 0
    # touching cache as a relevant header is newer
    FileUtils.touch cache_file
  end

  new_typemaps_deps = (buffer + typemaps_deps).uniq.sort
  if typemaps_deps != new_typemaps_deps
    File.open(cache_file, "w") { |f| f.puts new_typemaps_deps.join("\n")}

    # fake the modification time so that typemaps don't get regenerated on ALL fresh checkouts
    times = new_typemaps_deps.map { |h| File.mtime(h) }
    newest_time = times.max
    File.utime(newest_time, newest_time, cache_file)
  end

  exit(0)
end

if options[:scan]
  friends = `grep -r 'friend class Typemaps' #{xml_path}`.split("\n")
  friends.map! do |f|
    f =~ /^(.+class.+\.xml)\: /
    $1
  end
  friends.uniq!
  friends.map! do |f|
    File.open(f).read =~ /location file="(.*\.h)" line/
    $1
  end

  exit(0)
end

if options[:profiling]
  require 'ruby-prof'
  RubyProf.start
end

output_path = Pathname.new(ARGV[1] || "./")
header, source =
  TypemapGenerator.generate_forest(xml_path, basedir,
                                   options[:sloppy],
                                   options[:namespace],
                                   /#{options[:header_pattern]}/,
                                   options[:header_replacement],
                                   options[:macro_guard])
File.open(output_path + "typemaps.h",  "w").write(header)
File.open(output_path + "typemaps.#{options[:extension]}", "w").write(source)

if options[:profiling]
  profile = RubyProf.stop
  RubyProf::FlatPrinter.new(profile).print(STDOUT)
end
