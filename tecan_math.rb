module Tecan

	def Tecan.tipmask (aryOfTips) # takes an array of tips
		@tip_ar = aryOfTips
		x = 0
		@tip_ar.each do |a|
			x += 2**(a - 1)
		end
		return (x)
	end

	def Tecan.wellmask (aryOfWells, labwareWH) # takes an array of well positions and labware type
		@well_ar = aryOfWells
		@labwareWH = labwareWH
		@numWells = @labwareWH[0] * @labwareWH[1]

		# first 4 chr
		x = "%02X%02X" % [@labwareWH[0],@labwareWH[1]]

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


include Tecan
x = [1,2,3,4,5,6,7]
p(Tecan::wellmask(x, [12,8]))
