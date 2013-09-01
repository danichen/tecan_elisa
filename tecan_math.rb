module Tecan

	def tipmask (arg) # takes an array of tips
		@tip_ar = arg
		x = 0
		@tip_ar.each do |a|
			x += 2**(a - 1)
		end
		return (x)
	end

	def wellmask (arg1, arg2) # takes an array of well positions and labware type
		@well_ar = arg1
		@labware_ar = arg2
		@numWells = @labware_ar[0] * @labware_ar[1]

		# first 4 chr
		x = "%02X%02X" % [12,8]

		# the rest of the string
		
		(1..@numWells).each_slice(7) do |a|
			y = 48
			b = a & @well_ar
			b.each do |b|
				y = y + 2**a.index(b)
			end
		x = x + y.chr
		end
		return (x)
	end

end


#include Tecan
#x = [1, 5]
#p(Tecan::wellmask(x, [12,8]))
