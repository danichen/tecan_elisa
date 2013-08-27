#!/bin/env ruby
# encoding: utf-8

require "rubygems"
require "sequel"

# Global variables

DefaultVol = 60
MinVol = 10
MaxVol = 1000
MaxFactor = 100
StdConc = 200 # temp, will fetch it from database in future
Diluent = 'PTT' # temp, will fetch it from database in future


# command line arguments

argList = :migrate, :help

puts "Options: " << (ARGV.map!(&:to_sym)-argList).join(",") << " undefined and ignored"


# class definition

class DilutionSeries

  def initialize (series, sample)
    @series = series
    @sample = sample
    @numRep = @series.length/@series.uniq.length

		@dest_name = Destination.last.pk

    @sample_needed = 0
    @series.uniq.each {|a| @sample_needed += (DefaultVol/a)*@numRep}
    @sample_needed = 10 if  @sample_needed <= MinVol

    @output = Array.new
    self.resolve
    self.readout
		self.write_to_db
  end


  def resolve
    @series.uniq.sort.each_with_index do |x, i|
      i == 0 ? df = x : df = (x/@series.uniq.sort[i-1])
      df <= MaxFactor ? @output.push(x) : @output.push(factor(x,df))
    end
  end

  def factor (arg1, arg2)
    result = Array.new
    target = arg1
    df = arg2
    df_i = Hash.new
    
    (2..MaxFactor).each do |x|
      df_i.merge!(x => x - df/x) if df % x == 0
    end

    (df_i.select {|a,b| b > 0}).empty? ? best_df = MaxFactor : best_df = df_i.key(((df_i.select {|a,b| b > 0}).values.min))


    df/best_df <= MaxFactor ? result.push(target/df*best_df, target) : result.push(target/df*best_df, factor(target, df/best_df))

    return(result)
  end
 

  def readout
    p(@sample, @series.uniq.sort)
    p(@output.flatten)
    puts "Rep = #{@numRep}"
    p(@sample_needed.ceil)
    puts()
  end


  def write_to_db

		last_swell = Swell.last 	# start new scratch plate if new day or database empty
    if (last_swell.nil? or Date.parse(Destination.last.rundate) < Date.today)
			flag = 'even'
			well_index = 1
		else											# else just pick up where last dilution left off
			flag = last_swell.evenodd
			well_index = last_swell.position + 1
		end


		@output.flatten.each_with_index do |a, i|
			Swell.create do |b|
				b.evenodd = flag
				b.diluent = Diluent
				b.position = well_index
				b.destination_id = @dest_name
				if i == 0
					b.df = a
					b.dvolume = @sample_needed.ceil * (a - 1)
 					b.sample = @sample
					b.svolume = @sample_needed.ceil
				else
					b.df = a / @output.flatten[i-1]
					b.sample = Swell.last.position
					b.svolume = [(MaxVol / b.df).floor, @available_vol].min
					b.dvolume = b.svolume * (b.df - 1)
				end
			end


			# connect to dwells

			Dwell.where(:gname => @sample, :d1value => a).update(:swell_id => Swell.last.id)
			

			@available_vol = (Swell.last.dvolume + Swell.last.svolume) - DefaultVol * @numRep
			well_index += 1

  	end
  end

end



class Robot

	def initialize
		@dest_name = Destination.last.name
		p(@dest_name)
	end

	def dilution
		filename = @dest_name + ".dilution.gwl"
		f = File.open(filename, 'w')
			
		Swell.where(:destination_id => Destination.last.pk).each_with_index do |x, i|
			#f.puts [x.evenodd, x.diluent, x.dvolume, x.sample, x.svolume, x.position].join("\t")

			ci = i%8
			x.sample.is_a? Integer ? pos = x.sample : pos = 'BLAH'
			p(x.sample)

			f.puts ['A', 'Sample?', nil, nil, pos, nil, x.svolume, nil, nil, 2**(ci)].join(";")
			#f.puts ['D'
			#f.puts ['W'
		end


	f.close
	end
end

		


# connect to a database
DB = Sequel.connect('postgres://dchen:test1234@localhost/dilution')



# create tables

DB.drop_table(:destinations, :dwells, :swells, :cascade=>true)

DB.create_table? :destinations do
  primary_key :id
  String :name
  String :barcode
  String :rundate
  String :elisa
  Integer :event
  Integer :runseq
  Boolean :completed
  index :name
end

DB.create_table? :swells do
  primary_key :id, :on_delete => :cascade 
  String :evenodd
  String :diluent
	Float :df
  Integer :dvolume
  String :sample
  Integer :svolume
  Integer :position
	foreign_key :destination_id, :destinations
end

DB.create_table? :dwells do
  primary_key :id, :on_delete => :cascade 
  String :wlocation
  String :gname
  String :gtype
  String :sname
  String :d1name
  Float :d1value
  String :d1unit
  Integer :position
  Integer :volume
  foreign_key :destination_id, :destinations
	foreign_key :swell_id, :swells
end



# define relationship

class Destination < Sequel::Model
  one_to_many :dwells
	one_to_many :swells
end

class Dwell < Sequel::Model
  many_to_one :swells
  many_to_one :destinations
end

class Swell < Sequel::Model
  one_to_many :dwells
	many_to_one :destinations
end





# open and read softmax file


infile = File.open(ARGV[0].to_s, "rb:UTF-16LE")



Destination.create do |a| # spawn a new destination plate
  a.name = ARGV[0].to_s
  a.rundate = Date.today.to_s
end

manifest = Array.new # sample manifest

header = Hash.new # process file header

infile.gets.encode("UTF-8").split("\t").each_with_index do |a, i|
  header[a.gsub(/\s+/, "_").downcase.to_sym] = i.to_i
end

while line = infile.gets # read into database
  line = line.encode("UTF-8").split("\t")

  if !line[header[:descriptor1_value]].nil? then
    puts line.join("..")
    Dwell.create do |a|
      a.wlocation = line[header[:﻿well_location]] # weird encoding problem
      a.gname = line[header[:group_name]]
      a.gtype = line[header[:group_type]]
      a.sname = line[header[:sample_name]]
      a.d1name = line[header[:descriptor1_name]]
      a.d1value = line[header[:descriptor1_value]]
      a.d1unit = line[header[:descriptor1_units]]
      map = lambda{|x| x[0].ord-64+(x[1..-1].to_i-1)*8 }
      a.position = map.call(a.wlocation)
      a.volume = DefaultVol
      a.destination_id = Destination.last.pk
    end
    manifest.push(line[header[:group_name]]) # add sample to manifest
  end 
end
infile.close



# test area
Dwell.where(:gname => 'sample3').each do |a|
  p(a)
  p(Destination.where(:id => a.destination_id).first)
end

# make dilution series

manifest.uniq.each do |a|
  series = Array.new
  Dwell.where(:gname => a).each do |b|
    series.push(b.d1value)
  end
  dilutionSeries = DilutionSeries.new(series, a)
end

# output to file

dilution_script = Robot.new
dilution_script.dilution





Destination.each do |x|
  #puts x.rundate
  #puts x.values.values.to_a
end
