#!/bin/env ruby
# encoding: utf-8

# TO DO:	1)	handling concentration
#					2) 	beef up data model
#					3)	edge effect on even/odd
#					4)	large prime dilution factor (float df)
#					5)	reporting and error detection
#					6)	reorganize into more readable structure, informative var naming, etc.
#					7)	proper mixin of Tecan methods not calling it as module methods
#					8)	implement a "0" dilution factor

require "rubygems"
require "sequel"

require './tecan_math'

# Global variables

DefaultVol = 60
MinVol = 10
MaxVol = 1000
MaxFactor = 100
StdConc = 200 # temp, will fetch it from database in future (for handling concentration)
Diluent = 'PTT' # temp, will fetch it from database in future (for basic worklist)

$sample_list = Array.new


# command line arguments

argList = :migrate, :help

puts "Options: " << (ARGV.map!(&:to_sym)-argList).join(",") << " undefined and ignored"


# class definition

class DilutionSeries

  def initialize (series, sample)
    @series = series
    @sample = sample
    @numRep = @series.length/@series.uniq.length

		$sample_list.push(@sample)

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

		
			
		Swell.where(:destination_id => Destination.last.pk).each_slice(8) do |x|
			#puts [x.evenodd, x.diluent, x.dvolume, x.sample, x.svolume, x.position].join("\t")

			@series = Hash.new
			
			x.each_with_index do |y, i|
				@series.merge!({ y.id => (i+1)}) # might have issue with order?
			end
			
			paral = 0
			index = 1
			matrix = Array.new(8){Array.new(8)}

			x.each_with_index do |y, i|

				
				if (/^[0-9]*$/ === y.sample) # if dilution well
					matrix[paral][index] = { y.id => (i+1) }
				else
					paral += 1
					index = 0
					matrix[paral][index] = { y.id => (i+1) }
				end
				index += 1
			end

# WorkList Generation

			f.puts("#{AdvWL.pickup(1..@series.count)}\r\n")

			matrix.transpose.each_with_index do |x, i|

				# turn array of hash into hash
					collection = Hash.new
					x.compact.each do |y|
						p(y)
						collection.merge!(y)
					end
						
				
				#if first step, make a new worklist step, if not, only make new worklist step if there's
				# stuff to do

				if (i == 0)
					x.compact.each do |y|
						z = AdvWL.new(y)
						f.print("#{z.asps}\r\n")
					end
				
					whole = AdvWL.new(@series)
					f.print("#{whole.aspd}\r\n")
					y = AdvWL.new(collection)
					f.print("#{y.dis}\r\n")
					f.print("#{y.mix}\r\n")
				else
					if x.compact.length > 0
						y = AdvWL.new(collection)

						f.puts("#{y.asp}\r\n")
						f.puts("#{y.dis}\r\n")
						f.puts("#{y.mix}\r\n")
					end
	
				end
					
			end

			f.puts("#{AdvWL.drop}\r\n")

# (END) WorkList Generation

		end

	f.close
	end



# Pipet to final plate
	def final
		
		filename = @dest_name + ".final.gwl"
		f = File.open(filename, 'w')


		Dwell.order(:position).where(:destination_id => Destination.last.pk).each_slice(8) do |x|

			@tips = Array.new
			@src = Array.new
			@tgt = Array.new

			puts "=============================!"
			x.each_with_index do |y, i|
				
				if (y.d1value > 0)
					@tips.push(i+1)
					@tgt.push(y.position)
					@src.push(Swell.where(:id => y.swell_id).last.position)
					puts "#{y.position} : #{Swell.where(:id => y.swell_id).last.position}"
				end
			end

			p(@tips)
			p(@tgt)
			p(@src)



			# Begin Worklist Generation

			@target_pvol = Array.new(12,0) # standard sample pipet volume
			
			if (@tips.count > 0)
				@tips.each do |z|
					@target_pvol[z-1] = (DefaultVol - 10)
				end
				
				f.puts("#{AdvWL.pickup(@tips)}\r\n")

				#===========================Need consistent class method; can use all arrays other than initiate method (consider pulling that out of Command Writer)


				# Aspirate Well
				grid = ',23,2,1,'
				@loop = ',0,0);'
				liqCls = "Water free dispense for testing"

				f.puts ('B;Aspirate(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@src, [12,8])).to_s << @loop << "\r\n")



# Dispense
				grid = ',23,1,1,'

				f.puts ('B;Dispense(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@tgt, [12,8])).to_s << @loop << "\r\n")



				#===========================


				f.puts("#{AdvWL.drop}\r\n")
			

			end

		end
	end

end




class AdvWL

	include Tecan

	@@loop = ',0,0);'
	

	def initialize (lotsOfHash) # Do the collapsing of pipeting commands here
		@input = lotsOfHash
		print "new object input:"
		p(@input)



		@src_well_list = Array.new
		@tgt_well_list = Array.new
		@tip_list = Array.new
		@target_svol = Array.new(12,0)
		@target_dvol = Array.new(12,0)
		@target_tvol = Array.new(12,0)
		@target_mvol = Array.new(12,0)
		
		
		
		@input.each_key do |x| # Collapsing

			@record = Swell.where(:id => (x), :destination_id => Destination.last.pk)

			@target_svol[(@input[x] - 1)] += @record.last.svolume
			@target_dvol[(@input[x] - 1)] += @record.last.dvolume
			@target_tvol[(@input[x] - 1)] += (@record.last.svolume + @record.last.dvolume)
			@target_mvol[(@input[x] - 1)] += ((@record.last.svolume + @record.last.dvolume)*0.7).floor
			

			@tgt_well_list.push(@record.last.position)
			@src_well_list.push(@record.last.sample.to_i)
			@tip_list.push(@input[x])

		
		end

		p(@target_svol)
		p(@tip_list)
		p(src: @src_well_list)
		p(tgt: @tgt_well_list)

	end

# Aspirate Well
	def asp
		grid = ',23,2,1,'
		liqCls = "Water free dispense for testing"

		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@src_well_list, [12,8])).to_s << @@loop)
	end


# Aspirate Sample
  def asps
		grid = ',15,0,1,'
		liqCls = "Water free dispense for testing"

					
		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask([$sample_list.index(@record.last.sample)+1], [1,16])).to_s << @@loop)
	end


# Aspirate Diluent
	def aspd
		grid = ',22,0,1,'
		liqCls = "Water free dispense trough"

		return ('B;Aspirate(' << (Tecan.tipmask(1..@input.count)).to_s << ',"' << liqCls << '",' << @target_dvol.join(',') << grid << (Tecan.wellmask((1..@input.count).to_a, [1,8])).to_s << @@loop)
	end
		
# Dispense
	def dis
		grid = ',23,2,1,'
		liqCls = "Water free dispense for testing"

		return ('B;Dispense(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_tvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end

# Mix
	def mix
		times = ',1'
		grid = ',23,2,1,'
		liqCls = "Water free dispense for testing"
		return ('B;Mix(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_mvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << times << @@loop)
	end


# Pickup DiTi
	def AdvWL.pickup (arrTip)
		@tips = arrTip
		return ('B;GetDITI2(' << (Tecan.tipmask(@tips)).to_s << ',"DiTi 1000ul LiHa",0,0,0,70);')
	end

# Drop DiTi
	def AdvWL.drop
		return ('B;DropDITI(255,29,6,10,70,0);')
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
	Boolean :origin
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
#################### ADDING BACK EMPTY WELLS
	else
		Dwell.create do |a|
			a.wlocation = line[header[:﻿well_location]] # weird encoding problem
      a.gname = '0'
      a.gtype = '0'
      a.sname = '0'
      a.d1name = '0'
      a.d1value = '0'
      a.d1unit = '0'
      map = lambda{|x| x[0].ord-64+(x[1..-1].to_i-1)*8 }
      a.position = map.call(a.wlocation)
      a.volume = '0'
      a.destination_id = Destination.last.pk
		end
##################### END ADDING BACK
  end 
end
infile.close



# test area
Dwell.where(:gname => 'sample3').each do |a|
  p(a)
  p(Destination.where(:id => a.destination_id).first)
end

# make dilution series

manifest.uniq.each_with_index do |a, i|

	
  series = Array.new
  Dwell.where(:gname => a).each do |b|
    series.push(b.d1value)
  end
  dilutionSeries = DilutionSeries.new(series, a)
end

# output to file

dilution_script = Robot.new
dilution_script.dilution

puts "Sample\tPosition"
$sample_list.each_with_index do |x, i|
	puts "#{x}\t#{i+1}"
end

# make final plate

dilution_script.final



