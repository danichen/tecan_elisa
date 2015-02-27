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

Final_Pipet_Vol = 100

DefaultVol = 150 #(default = 120) : Minimum volume in the S-Wells for each replicate (use 250 to make spiking easy)
PreSpike_DefaultVol =200
PostSpike_DefaultVol =200

Non_Series_Vol = (Final_Pipet_Vol + 30) * 2 # Volume to use in non-series S-wells for each replicate (default= + 20)

MinVol = 20 # Minimum Pipeting Volume, suggest 20 or more for Tecan (use 16 for 800 MaxVol)

Cutoff_Vol = 150 # Volume below which a change to p200 is required

Mix_Rep = 2

Large_Mixing_Vol = 300


MaxVol = 800 # scratch well capacity, default = 1000, trimmed to 500/800 improve liquid detection


MaxFactor = (MaxVol / MinVol).floor

StdConc = 200 # temp, will fetch it from database in future (for handling concentration)
Diluent = 'PTT' # temp, will fetch it from database in future (for basic worklist)

$sample_list = Array.new

$spike_hash = {'None' => 13, 'Low' => 14, 'Mid' => 15, 'Hi' => 16}


# Change scratch well starting position
puts "Enter column index on scratch plate to start the dilutions:"
@@starting_column = STDIN.gets.chomp.to_i
if (@@starting_column < 2 || @@starting_column > 12)
	@@starting_column = 1
end




# command line arguments

argList = :migrate, :help

puts "Options: " << (ARGV.map!(&:to_sym)-argList).join(",") << " undefined and ignored"



# class definition

class DilutionSeries

  def initialize (series, sample, spike)

		@spike_or_not = spike

    @series = series
    @sample = sample
    @numRep = @series.length/@series.uniq.length

		@new_sample_flag = 1

		@dest_name = Destination.last.pk

    @output = Array.new
    self.resolve
    self.readout
		self.write_to_db
  end


  def resolve
    @series.uniq.sort.each_with_index do |x, i|
      if (i == 0)

				df = x
				@correction_factor = 1

				Swell.all.each do |s|   # Check if base_sample is already in use
					if s.base_sample == $base_sample_hash[@sample] && s.dilution == x
						@new_sample_flag = 0
						puts("not a new sample do not factor")
						@correction_factor = df/s.df
						puts("#{@correction_factor}")
						break
					end
				end
				
				if (@new_sample_flag != 0)
					df <= MaxFactor ? @output.push(x) : @output.push(factor(x,df))
				else
					@output.push(x)
				end

			else
				df = (x/@series.uniq.sort[i-1])
      	df <= MaxFactor ? @output.push(x) : @output.push(factor(x,df))
			end
    end

		use_vol = DefaultVol
		use_vol = PostSpike_DefaultVol if @spike_or_not == 2
		use_vol = PreSpike_DefaultVol if @spike_or_not == 1

		@sample_needed = 0
    @output.flatten.uniq.each {|a| @sample_needed += (use_vol/(a))*@correction_factor*@numRep}
    @sample_needed = MinVol if  @sample_needed <= MinVol

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
			well_index = (@@starting_column - 1) * 8 + 1
		else											# else just pick up where last dilution left off
			flag = last_swell.evenodd
			well_index = last_swell.position + 1
		end


		@output.flatten.each_with_index do |a, i|
			Swell.create do |b|
				if (a == @series.uniq.sort.first)
					b.first = 'y'
					b.prespike = Dwell.where(:gname => @sample, :working_d1value => a).last.prespike
				end

				if Dwell.where(:gname => @sample, :working_d1value => a).last != nil
					b.postspike = Dwell.where(:gname => @sample, :working_d1value => a).last.postspike
				end

				b.evenodd = flag
				b.position = well_index
				b.destination_id = @dest_name
				b.base_sample = $base_sample_hash[@sample]

				b.dilution = a

				if i == 0
					
					if (@new_sample_flag != 0) # if first well is a new dilution

						b.df = a
						b.dvolume = @sample_needed.ceil * (a - 1)
 						b.sample = @sample
						b.svolume = @sample_needed.ceil
					
						### Use fixed volume for non-series S-wells
						b.svolume = Non_Series_Vol if @output.flatten.length == 1 && b.df == 1 && b.prespike == nil && b.postspike == nil

						b.vol_left = b.dvolume + b.svolume ###(new)
						b.tvolume = b.vol_left

					else # if first well is taken from well from a previous dilution series

						######### Find the matching dilution and populate b.* #############

						Swell.all.each do |s|   # Check if base_sample is already in use
							if (s.base_sample == b.base_sample && s.dilution == b.dilution)
								b.df = s.df
								#b.dvolume = s.dvolume
								b.dvolume = @sample_needed.ceil * (b.df - 1)
								b.sample = s.sample
								#b.svolume = s.svolume
								b.svolume = @sample_needed.ceil

								### Use fixed volume for non-series S-wells
								b.svolume = Non_Series_Vol if @output.flatten.length == 1 && b.df == 1 && b.prespike == nil && b.postspike == nil
								
								b.vol_left = b.dvolume + b.svolume
								b.tvolume = b.vol_left

								if (/^[0-9]*$/ === b.sample) # if sample is not from real sample source
									vol_left = Swell.where(:position => b.sample).first.vol_left - b.svolume
									Swell.where(:position => b.sample).first.update(:vol_left => vol_left)
								end
								
								break
							end
						end
						
					end

				else
					b.df = a / @output.flatten[i-1]
					b.sample = Swell.last.position
					b.svolume = [(MaxVol / b.df).floor, @available_vol].min
					b.dvolume = b.svolume * (b.df - 1)
					b.vol_left = b.dvolume + b.svolume
					b.tvolume = b.vol_left

					vol_left = Swell.last.vol_left - b.svolume
					Swell.last.update(:vol_left => vol_left)

				end
			end


			# connect to dwells

			Dwell.where(:gname => @sample, :working_d1value => a).update(:swell_id => Swell.last.id)
			

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

		
			
		Swell.where(:destination_id => Destination.last.pk).order(:id).each_slice(8) do |x|

			@series = Hash.new
			@series_s200 = Hash.new

			@series_s1000 = Hash.new
			@series_diluent = Hash.new

			diluent_required = 0
			
			x.each_with_index do |y, i|
				@series.merge!({ y.id => (i+1)}) # might have issue with order?

				if y.dvolume != 0
					@series_s200.merge!({ y.id => (i+1)}) if y.svolume < Cutoff_Vol
					@series_s1000.merge!({ y.id => (i+1)}) if y.svolume >= Cutoff_Vol
					@series_diluent.merge!({ y.id => (i+1)})
				else
					@series_s200.merge!({ y.id => (i+1)}) if y.svolume < Cutoff_Vol
					@series_s1000.merge!({ y.id => (i+1)}) if y.svolume >= Cutoff_Vol
				end
				
				diluent_required = diluent_required + y.dvolume
			end

			puts ("Need diluent #{diluent_required} uL")
			
			paral = 0
			index = 1
			matrix = Array.new(8){Array.new(8)}

			x.each_with_index do |y, i|

				
				if (/^[0-9]*$/ === y.sample) # if dilution well
					
					(index-1).times do |z|
						if matrix[paral][(z+1)] != nil && Swell.where(:id => matrix[paral][(z+1)].keys).first.sample == y.sample

							index = z+1
							paral += 1
							break
						end
					end

					matrix[paral][index] = { y.id => (i+1) }
				else
					index = 0
					matrix[paral][index] = { y.id => (i+1) }
					paral += 1
				end
				index += 1
			end


# WorkList Generation (step-wise in a series contained within a column of 8 wells)

			whole = AdvWL.new(@series)
			whole_s200 = AdvWL.new(@series_s200)
			whole_s1000 = AdvWL.new(@series_s1000)
			whole_diluent = AdvWL.new(@series_diluent)

			f.puts("C;===Pipetting wells #{@series.keys} with tips #{@series.values}===\r\n")

			if (diluent_required > 0)

				f.puts("#{whole_diluent.pickupspecific_1000}\r\n") # pick up p1000 for all wells needing diluent
				f.print("#{whole_diluent.asptrough}\r\n")  # aspirate required diluent for all eligile wells
				f.print("#{whole_diluent.disd}\r\n") # dispense diluent into S-wells
				f.print("#{whole_diluent.dropspecific}\r\n") # drop all diluent tips
			end

			f.puts("#{whole_s200.pickupspecific_200}\r\n")
			f.puts("#{whole_s1000.pickupspecific_1000}\r\n")

			matrix.transpose.each_with_index do |x, i| # Turn series into parallel steps

				# flatten an array of hash into a hash   #### <<<MAKE SOME TIP DECISION HERE>>> ####
					collection = Hash.new
					
					need_tip = Hash.new #(new)
					need_small = Hash.new #(new)
					need_large = Hash.new #(new)
					need_large_mixing = Hash.new #(new)

					need_wash = Hash.new #(new)

					need_prespike = Hash.new #(new/new)
					need_prespike_then_large = Hash.new #(new/new)
					
					x.compact.each do |y|
		
						collection.merge!(y)


						if Swell.where(:id => y.keys).last.svolume > 0

							need_wash.merge!(y) if Swell.where(:id => y.keys).last.dvolume > 0
				
							if Swell.where(:id => y.keys).last.svolume < Cutoff_Vol
								puts "have to use p200"
								need_small.merge!(y)
								need_large_mixing.merge!(y) if Swell.where(:id => y.keys).last.tvolume > Large_Mixing_Vol
							else
								puts "stick with p1000"
								need_large.merge!(y)
								need_large_mixing.merge!(y) if Swell.where(:id => y.keys).last.prespike != nil && Swell.where(:id => y.keys).last.tvolume > Large_Mixing_Vol
							end
							need_tip.merge!(y)
							
						end

						if Swell.where(:id => y.keys).last.prespike != nil
							puts "need prespike"
							need_prespike.merge!(y)
						end
					end
						
				
				#if first step, make a new worklist step, if not, only make new worklist step if there's
				# stuff to do

				# 'whole' (declared above outside the steps) is every well in a column(8) that has a volume (from @series, a hash)
				# 'step' is every well in a column(8) that can be pipetted at the same step
				# 'small' is every well in 'step' that needs a tip change to p200 to pipet from sample source
				# 'large' is every well in 'step' that needs a tip change to p1000 to pipet from sample source
				# 'both' is the union of 'small' and 'large'
				#		>> small and large don't include pipettor where d-volume = 0, because s-volume then must be larger than cutoff
				# 	>> and the tips aren't contaminated

				step = AdvWL.new(collection)
				small = AdvWL.new(need_small)
				large = AdvWL.new(need_large)
				both = AdvWL.new(need_tip)

				pre_spike = AdvWL.new(need_prespike)

				large_mixing = AdvWL.new(need_large_mixing)

				wash = AdvWL.new(need_wash)

				

				if (i == 0)  # i.e., pipetter will aspirate from sample source

				
					if need_tip.length > 0 # Possible bug here?

						f.puts("#{both.asps_compact_smart_2}")
						
						f.print("#{both.diss}\r\n")
						f.print("#{wash.mix_200}\r\n") if (need_wash.length > 0)
						f.puts("#{both.dropspecific}\r\n")
					end


					### DO PRESPIKING HERE:
					if (need_prespike.length > 0)
						f.puts("C;Pre-spiking===\r\n")
						f.puts("#{pre_spike.pickupspecific_200}\r\n")

						# Aspirate spike, unpacked in case multiple spike series occur in the same octet
						need_prespike.each do |y|
							z = AdvWL.new({y[0] => y[1]}) # Hash turned into arrays of array by '.each'
							f.puts("#{z.asp_spike}\r\n")
						end
						# end Aspirate spike

						f.print("#{pre_spike.dis_spike}\r\n") # packed dis_spike
						f.print("#{pre_spike.mix_spike}\r\n") # packed mix_spike
						f.puts("#{pre_spike.dropspecific}\r\n")
						f.puts("C;End pre-spiking===\r\n")
					end

					
					# Swap out p200 for p1000
					
					if need_large_mixing.length > 0
						f.puts("#{large_mixing.pickupspecific_1000}\r\n")
						f.print("#{large_mixing.mix_1000}\r\n")
						f.puts("#{large_mixing.dropspecific}\r\n")
					end


				else # Pipet from a S-well
					if x.compact.length > 0
						

						f.puts("#{both.asp_smart}") # Need to unpack if multi-aspirating from a single sample
						f.print("#{both.diss}\r\n")
						f.puts("#{both.mix_200}\r\n")
						f.puts("#{both.dropspecific}\r\n")
		

						### DO PRESPIKING HERE:
						if (need_prespike.length > 0)
							f.puts("C;Pre-spiking===\r\n")
							f.puts("#{pre_spike.pickupspecific_200}\r\n")

							# Aspirate spike, unpacked in case multiple spike series occur in the same octet
							need_prespike.each do |y|
								z = AdvWL.new({y[0] => y[1]}) # Hash turned into arrays of array by '.each'
								f.puts("#{z.asp_spike}\r\n")
							end
							# end Aspirate spike

							f.print("#{pre_spike.dis_spike}\r\n") # packed dis_spike
							f.print("#{pre_spike.mix_spike}\r\n") # packed mix_spike
							f.puts("#{pre_spike.dropspecific}\r\n")
							f.puts("C;End pre-spiking===\r\n")
						end


						if need_large_mixing.length > 0
							f.puts("#{large_mixing.pickupspecific_1000}\r\n")
							f.print("#{large_mixing.mix_1000}\r\n")
							f.puts("#{large_mixing.dropspecific}\r\n")
						end

					end
	
				end
					
			end


		# (END) WorkList Generation

		end

		# DO POST-SPIKE HERE

		post_spike = Array.new

		Swell.where(:destination_id => Destination.last.pk).order(:id).each do |x|

			if x.postspike != nil # SOFTMAX: change to postspike (check!)
				post_spike.push(x)
			end
		end

		f.puts("C;Do post spiking.\r\n")

		post_spike.each_slice(8) do |x|

			post_tips = Array.new
			post_svol = Array.new(12,0)
			post_mvol = Array.new(12,0)
			post_tgt = Array.new

			x.length.times do |y|
				post_tips.push(y+1)
			end

			if x.length > 0
				f.puts("#{AdvWL.pickup200(post_tips)}\r\n")

				post_liqCls = "Scratch"
				post_loop = ',0,0);'
				post_grid = ',13,0,1,'

				x.each_with_index do |y, i|

					w = y.vol_left * 0.02
					w = w.to_i if w % 1 == 0

					svol = Array.new(12,0)
					svol[i] = w

					post_svol[i] = w ### COMFIRM
					post_mvol[i] = 160


					post_tgt.push(y.position)

					tgt = Array.new
					tgt.push($spike_hash[y.postspike]) # SOFTMAX: change to postspike (check!)
					

					# Aspirate post-spike one-by-one

					f.puts('B;Aspirate(' << (Tecan.tipmask(post_tips[i,1])).to_s << ',"' << post_liqCls << '",'<< svol.join(',') << post_grid << (Tecan.wellmask(tgt, [1,16])).to_s << post_loop << "\r\n")
					

				end

				# Dispense post-spike together

				post_liqCls = "Scratch"
				post_grid = ',16,1,1,'

				f.puts('B;Dispense(' << (Tecan.tipmask(post_tips)).to_s << ',"' << post_liqCls << '",'<< post_svol.join(',') << post_grid << (Tecan.wellmask(post_tgt, [12,8])).to_s << post_loop << "\r\n")

				# Mix 3 times
				
				post_liqCls = "FastMixing"
				post_times = ",3"

				f.puts('B;Mix(' << (Tecan.tipmask(post_tips)).to_s << ',"' << post_liqCls << '",'<< post_mvol.join(',') << post_grid << (Tecan.wellmask(post_tgt, [12,8])).to_s << post_times << post_loop << "\r\n")

				# Dump tips

				f.puts("#{AdvWL.drop}\r\n")
				
				
			end
		end



	f.close
	end



# Pipet to final plate
	def final

		@last_src = Array.new
		
		filename = @dest_name + ".final.gwl"
		f = File.open(filename, 'w')

		dispense_buffer = ""
		aspirate_buffer = ""
		diti_buffer = ""
		rep = 1

		mix_buffer = "" # new


		Dwell.order(:position).where(:destination_id => Destination.last.pk).each_slice(8) do |x|

			@tips = Array.new
			@src = Array.new
			@tgt = Array.new

			@mix = Array.new

			puts "\n=============================!"
			x.each_with_index do |y, i|
				
				
				if (y.d1value > 0) # y.working.d1value
					@tips.push(i+1)
					@tgt.push(y.position)
					@src.push(Swell.where(:id => y.swell_id).last.position)
					puts "#{y.position} : #{Swell.where(:id => y.swell_id).last.position}"

				end
			end

			@target_pvol = Array.new(12,0) # standard sample pipet volume
			@mix_vol = Array.new(12,0)

			@condition_vol = Array.new(12,0)


			
			if @src == @last_src # if just a replicate:

				puts "duplicate!"
				rep += 1; p(rep)
				liqCls = "Sample Transfer"
				@loop = ',0,0);'

				@tips.each do |z|
					@target_pvol[z-1] = Final_Pipet_Vol
				end

				# Dispense
				grid = ',16,2,1,'
				dispense_buffer = dispense_buffer << ('B;Dispense(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@tgt, [12,8])).to_s << @loop << "\r\n") if @tgt.length > 0

				# Aspirate
				grid = ',16,1,1,'

				@tips.each do |z|
					@target_pvol[z-1] = (Final_Pipet_Vol * rep) + 20
				end
			

				aspirate_buffer = ('B;Aspirate(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@src, [12,8])).to_s << @loop << "\r\n")

			### Test to improve duplicate accuracy

				@tips.each do |z|
					@condition_vol[z-1] = 10
				end
				aspirate_buffer = aspirate_buffer << ('B;Dispense(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @condition_vol.join(',') << grid << (Tecan.wellmask(@src, [12,8])).to_s << @loop << "\r\n")
			###

				

			else # else, current column is not another replicate
				puts "new replicate set: write previous!"

				

				if dispense_buffer != ""
					f.print(diti_buffer)
					f.print(aspirate_buffer)
					f.print(dispense_buffer)
					f.puts("#{AdvWL.drop}\r\n")
				end

				rep = 1 #initialize
				aspirate_buffer = ""
				dispense_buffer = ""
				diti_buffer = ""
				
				mix_buffer = ""

				# MIX (new)

				@tips.each_with_index do |z, i|
					p(@src)
					temp_vol = 0
					if (Swell.where(:position => @src[i]).last.vol_left * 0.5).floor.to_i > 180
						temp_vol = 180
					else
						temp_vol = (Swell.where(:position => @src[i]).last.vol_left * 0.5).floor.to_i
					end

					# BSA uses p1000 to distribute diluted sample
					temp_vol = (Swell.where(:position => @src[i]).last.vol_left * 0.5).floor.to_i
					
					
					@mix_vol[z-1] = temp_vol if z != 0
				end
				
				times = ",#{Mix_Rep}"
				grid = ',16,1,1,'
				liqCls = "FastMixing"
				@loop = ',0,0);'

	
				mix_buffer = ('B;Mix(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @mix_vol.join(',') << grid << (Tecan.wellmask(@src, [12,8])).to_s << times << @loop << "\r\n") 
	

				# Dispense

				@tips.each do |z|
					@target_pvol[z-1] = Final_Pipet_Vol
				end


				liqCls = "Sample Transfer"
				grid = ',16,2,1,'

				
				dispense_buffer = dispense_buffer << ('B;Dispense(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@tgt, [12,8])).to_s << @loop << "\r\n") if @tgt.length > 0

				# Aspirate
				grid = ',16,1,1,'

				@tips.each do |z|
					@target_pvol[z-1] = (DefaultVol * rep)
				end
			

				aspirate_buffer = ('B;Aspirate(' << (Tecan.tipmask(@tips)).to_s << ',"' << liqCls << '",'<< @target_pvol.join(',') << grid << (Tecan.wellmask(@src, [12,8])).to_s << @loop << "\r\n")

				# Diti

				diti_buffer = ("#{AdvWL.pickup(@tips)}\r\n")


				p(@src)
				p(@last_src)
				@last_src = Array.new(@src)
					
				#===========================Need consistent class method; can use all arrays other than initiate method (consider pulling that out of Command Writer)
			end
			
		end

		###### Print last buffered column...NOT DRY BIG TIME...:(

		if dispense_buffer != ""
			f.print(diti_buffer)
			f.print(aspirate_buffer)
			f.print(dispense_buffer)
			f.puts("#{AdvWL.drop}\r\n")
		end

	f.close

	end



end




class AdvWL

	include Tecan

	@@loop = ',0,0);'
	

	def initialize (lotsOfHash) # Do the collapsing of pipeting commands here
		@input = lotsOfHash

		if @input != nil

		end

		@src_sample_list = Array.new #(new)

		@src_well_list = Array.new
		@tgt_well_list = Array.new
		@tip_list = Array.new
		@target_svol = Array.new(12,0)
		@target_dvol = Array.new(12,0)
		@target_tvol = Array.new(12,0)
		@target_mvol = Array.new(12,0)
		@target_mvol_mixtip = Array.new(12,0)
		
		
		
		@input.each_key do |x| # Collapsing

			@record = Swell.where(:id => x, :destination_id => Destination.last.pk)

			@target_svol[(@input[x] - 1)] += @record.last.svolume
			@target_dvol[(@input[x] - 1)] += @record.last.dvolume
			@target_tvol[(@input[x] - 1)] += (@record.last.svolume + @record.last.dvolume)

			# Mixing volume used if p200 are mixed in, else use 'tvol'
			@record.last.svolume <= Cutoff_Vol ? temp_vol = 200 : temp_vol = (@record.last.svolume + @record.last.dvolume)
			@target_mvol_mixtip[(@input[x] - 1)] += (temp_vol * 0.8).floor

			@target_mvol[(@input[x] - 1)] += ((@record.last.svolume + @record.last.dvolume) * 0.8).floor
			

			@tgt_well_list.push(@record.last.position)
			@src_well_list.push(@record.last.sample.to_i)
			@tip_list.push(@input[x])

			@src_sample_list.push($base_sample_hash.values.compact.uniq.index($base_sample_hash[@record.last.sample]).to_i+1) #new/new
			

		
		end

		if @input.count != 0
			print "\n============\nnew object input(well_pos => tip):"
			#p(caller_locations(1,1)[0].label)
			#p(caller[0][/`([^']*)'/, 1])
			p(@input)
			p(@target_svol)
			p(@tip_list)
			p(src: @src_well_list)
			p(tgt: @tgt_well_list)
		end

	end


# Aspirate Well
	def asp
		grid = ',16,1,1,'
		liqCls = "Scratch"

		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@src_well_list, [12,8])).to_s << @@loop)
	end


# Aspirate Well together, but unpack if multiple tips are drawing from the same source
	def asp_smart
		grid = ',16,1,1,'
		liqCls = "Scratch"

		return_stack = ""

		if @tip_list.length != @src_well_list.uniq.length


			@tip_list.each_with_index do |x, i|

				svol_array = Array.new(12,0)
				svol_array[(x-1)] = @target_svol[(x-1)]
				
				tip_list = Array.new
				tip_list[0] = @tip_list[i]

				well_list = Array.new
				well_list[0] = @src_well_list[i]
				

				return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list)).to_s << ',"' << liqCls << '",'<< svol_array.join(',') << grid << (Tecan.wellmask(well_list, [12,8])).to_s << @@loop << "\r\n")
			end
		else
			return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@src_well_list, [12,8])).to_s << @@loop << "\r\n")
		end

		return (return_stack)

	end


# Aspirate Sample
  def asps
		grid = ',14,0,1,'
		liqCls = "Sample Transfer"
					
		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask([$sample_list.compact.uniq.index(@record.last.sample)+1], [1,16])).to_s << @@loop)
	end

# Aspirate Sample together (new)
  def asps_compact
		grid = ',14,0,1,'
		liqCls = "Sample Transfer"
					
		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@src_sample_list, [1,16])).to_s << @@loop)
	end



# Aspirate Sample together but not if a bunch of tips are drawing from the same sample (new)
  def asps_compact_smart
		grid = ',14,0,1,'
		liqCls = "Sample Transfer"

		return_stack = ""

		# if @tip_list.length != @src_sample_list.uniq.length
		if (1)



			@tip_list.each_with_index do |x, i|

				svol_array = Array.new(12,0)
				svol_array[(x-1)] = @target_svol[(x-1)]
				
				tip_list = Array.new
				tip_list[0] = @tip_list[i]

				well_list = Array.new
				well_list[0] = @src_sample_list[i]
				

				return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list)).to_s << ',"' << liqCls << '",'<< svol_array.join(',') << grid << (Tecan.wellmask(well_list, [1,16])).to_s << @@loop << "\r\n")
			end
		else
			return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@src_sample_list, [1,16])).to_s << @@loop << "\r\n")
		end



		return (return_stack)
	end



# Aspirate Sample together but not if a bunch of tips are drawing from the same sample, accomodates sample number > 16 (new)
  def asps_compact_smart_2
		grid1 = ',14,0,1,'
		grid2 = ',15,0,1,'
		liqCls = "Sample Transfer"

		return_stack = ""

		if @tip_list.length != @src_sample_list.uniq.length # More tip than samples, break into single pipeting


			@tip_list.each_with_index do |x, i|

				svol_array = Array.new(12,0)
				svol_array[(x-1)] = @target_svol[(x-1)]
				
				tip_list = Array.new
				tip_list[0] = @tip_list[i]

				well_list = Array.new
				if @src_sample_list[i] <= 16
					well_list[0] = @src_sample_list[i]
				
					return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list)).to_s << ',"' << liqCls << '",'<< svol_array.join(',') << grid1 << (Tecan.wellmask(well_list, [1,16])).to_s << @@loop << "\r\n")
				else # Pipet from second sample carrier
					well_list[0] = @src_sample_list[i]-16

					return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list)).to_s << ',"' << liqCls << '",'<< svol_array.join(',') << grid2 << (Tecan.wellmask(well_list, [1,16])).to_s << @@loop << "\r\n")
				end
				
			end
		else # Pipet together, but break sample list into sample carrier 1 and sample carrier 2
			svol_array1 = Array.new(12,0)
			svol_array2 = Array.new(12,0)
			tip_list1 = Array.new
			tip_list2 = Array.new
			well_list1 = Array.new
			well_list2 = Array.new
			
			@tip_list.each_with_index do |x, i|
				if @src_sample_list[i] <= 16
					svol_array1[(x-1)] = @target_svol[(x-1)]
					tip_list1.push(@tip_list[i])
					well_list1.push(@src_sample_list[i])
				else
					svol_array2[(x-1)] = @target_svol[(x-1)]
					tip_list2.push(@tip_list[i])
					well_list2.push(@src_sample_list[i]-16)
				end
			end

				if tip_list1.count > 0
					return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list1)).to_s << ',"' << liqCls << '",'<< svol_array1.join(',') << grid1 << (Tecan.wellmask(well_list1, [1,16])).to_s << @@loop << "\r\n")
				end

				if tip_list2.count > 0
					return_stack = return_stack << ('B;Aspirate(' << (Tecan.tipmask(tip_list2)).to_s << ',"' << liqCls << '",'<< svol_array2.join(',') << grid2 << (Tecan.wellmask(well_list2, [1,16])).to_s << @@loop << "\r\n")
				end
			
		
		end



		return (return_stack)
	end



# Aspirate Diluent from Trough
	def asptrough
		grid = ',22,4,1,'
		liqCls = "OverScratch"

		
		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",' << @target_dvol.join(',') << grid << (Tecan.wellmask(@tip_list, [1,8])).to_s << @@loop)
	end

# Aspirate diluent volume from scratch well (new)
	def aspd
		grid = ',16,1,1,'
		liqCls = "Scratch"

		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_dvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end
		
# Dispense total volume (change)
	def dist
		grid = ',16,1,1,'
		liqCls = "Scratch"

		return ('B;Dispense(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_tvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end

# Dispense total volume (new)
	def diss
		grid = ',16,1,1,'
		liqCls = "Scratch"

		return ('B;Dispense(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_svol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end

# Aspirate spike sample (new)##############

	def asp_spike
		#grid = ',14,0,1,'
		grid = ',13,0,1,'
		liqCls = "Scratch"

		spike_vol = Array.new
		@target_tvol.each do |v|
			w = v * 0.02
			w = w.to_i if w % 1 == 0
			spike_vol.push(w)
		end

		coord = Array.new
		coord.push($spike_hash[Swell.where(:position => @tgt_well_list[0]).last.prespike])
					
		return ('B;Aspirate(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< spike_vol.join(',') << grid << (Tecan.wellmask(coord, [1,16])).to_s << @@loop)
	end


# Dispense spike volume (new)
	def dis_spike
		grid = ',16,1,1,'
		liqCls = "Scratch"

		spike_vol = Array.new
		@target_tvol.each do |v|
			w = v * 0.02
			w = w.to_i if w % 1 == 0
			spike_vol.push(w)
		end

		return ('B;Dispense(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< spike_vol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end

# Mix spike volume
	def mix_spike
		times = ",#{Mix_Rep}"
		grid = ',16,1,1,'
		liqCls = "FastMixing"

		spike_mix_vol = Array.new
		@target_tvol.each do |v|
			(v < 100) ? w = v : w = 100
			spike_mix_vol.push(w)
		end
		return ('B;Mix(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< spike_mix_vol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << times << @@loop)
	end


# Dispense diluent volume in sratch plate (new)
	def disd
		grid = ',16,1,1,'
		liqCls = "OverScratch"

		return ('B;Dispense(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_dvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << @@loop)
	end

# Mix with set of only p1000 tips
	def mix_1000
		times = ",#{Mix_Rep}"
		grid = ',16,1,1,'
		liqCls = "FastMixing"
		return ('B;Mix(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_mvol.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << times << @@loop)
	end

# Mix with p200 included in the tip set
	def mix_200
		times = ",#{Mix_Rep}"
		grid = ',16,1,1,'
		liqCls = "FastMixing"
		return ('B;Mix(' << (Tecan.tipmask(@tip_list)).to_s << ',"' << liqCls << '",'<< @target_mvol_mixtip.join(',') << grid << (Tecan.wellmask(@tgt_well_list, [12,8])).to_s << times << @@loop)
	end


# Pickup DiTi
	def AdvWL.pickup (arrTip)
		@tips = arrTip
		return ('B;GetDITI2(' << (Tecan.tipmask(@tips)).to_s << ',"DiTi 1000ul SBS LiHa",0,0,0,70);')
	end

# Pickup DiTi
	def AdvWL.pickup200 (arrTip)
		@tips = arrTip
		return ('B;GetDITI2(' << (Tecan.tipmask(@tips)).to_s << ',"DiTi 200ul SBS LiHa",0,0,0,70);')
	end

# Drop Specific DiTi (new)
	def dropspecific
		return ('B;DropDITI(' << (Tecan.tipmask(@tip_list)).to_s <<  ',22,6,10,70,0);')
	end

# Pickup Specific p1000 DiTi
	def pickupspecific_1000
		return ('B;GetDITI2(' << (Tecan.tipmask(@tip_list)).to_s << ',"DiTi 1000ul SBS LiHa",0,0,0,70);')
	end

# Pickup Specific p200 DiTi
	def pickupspecific_200
		return ('B;GetDITI2(' << (Tecan.tipmask(@tip_list)).to_s << ',"DiTi 200ul SBS LiHa",0,0,0,70);')
	end

# Drop All DiTi
	def AdvWL.drop
		return ('B;DropDITI(255,22,6,10,70,0);')
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
	Float :df
  Integer :dvolume
  String :sample
  Integer :svolume
  Integer :position
	Boolean :origin
	Integer :vol_left
	Integer :tvolume
	String :first
	String :base_sample
	Float :dilution
	String :prespike
	String :postspike
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
	Float :working_d1value #holds working dilution factors 
  String :d1unit
  Integer :position
  Integer :volume
	String :base_sample
	String :prespike
	String :postspike
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
			temp = a.gname.split("-")
			a.base_sample = temp[0] #reverse for the new SOFTMAX format sample:spike (check!)
			a.prespike = temp[1]
			temp = a.gname.split(":")
			a.postspike = temp[1]
    end

#################### ADDING BACK EMPTY WELLS
	else
		Dwell.create do |a|
			a.wlocation = line[header[:﻿well_location]] # weird encoding problem
      a.gname = nil
      a.gtype = '0'
      a.sname = '0'
      a.d1name = '0'
      a.d1value = '0'
      a.d1unit = '0'
      map = lambda{|x| x[0].ord-64+(x[1..-1].to_i-1)*8 }
      a.position = map.call(a.wlocation)
      a.volume = '0'
      a.destination_id = Destination.last.pk
			a.base_sample = nil
			a.prespike = nil
			a.postspike = nil
		end
##################### END ADDING BACK
  end 
end
infile.close


# make dilution series, Swell ordered by Dwell position to optimize S->D transfer

Dwell.order(:position).where(:destination_id => Destination.last.pk).each do |a|
	$sample_list.push(a.gname)
end


# Make a hash to translate variant sample name to base_sample name

$base_sample_hash = Hash.new

$sample_list.compact.uniq.each do |a|
	b = (a.split("-")[0]).split(":")[0] ################ CHANGE when SOFTMAX updated (check!)
	$base_sample_hash.merge!(a => b)
end

puts "=================="
p($base_sample_hash)
puts "=================="


# Process each sample variants here

$sample_list.compact.uniq.each do |a|
	
  series = Array.new
	conc = 0

	# Get stock concentration from user
	if (p(Dwell.where(:gname => a).last.d1name) == 'Concentration')
		puts "Enter stock concentration for #{a} (ng/mL):"
		conc = STDIN.gets.chomp.to_f
	end

	# Pass series value to DilutionSeries class; convert conc to df if necessary
  Dwell.where(:gname => a).each do |b|
		if conc > 0
			b.working_d1value = conc/b.d1value
		else
			b.working_d1value = b.d1value
		end
		Dwell.where(:id => b.id).update(:working_d1value => b.working_d1value)
    series.push(b.working_d1value)
  end

	# Determine whether the series is spiked
	spike = 0

	Dwell.where(:gname => a).each do |b|
		if b.prespike != nil || b.postspike != nil
			spike = 1 if b.prespike != nil
			spike = 2 if b.postspike != nil
			break
		end
	end

  dilutionSeries = DilutionSeries.new(series, a, spike)
end

# output to file

dilution_script = Robot.new
dilution_script.dilution


# output sample location assignment
puts "Sample\tPosition"
$sample_list.compact.uniq.each_with_index do |x, i|
	puts "#{x}\t#{i+1}"
end

puts "================" # Print unique samples

$base_sample_hash.values.compact.uniq.each_with_index do |x, i|
	puts "#{x}\t#{i+1}"
end

# make final plate

dilution_script.final



