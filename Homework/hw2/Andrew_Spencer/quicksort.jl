# QUICKSORT FOR A traditional row-based array with multiple columns
function quicksort2D(a)#,left,right)
    
    lena = length(a[:,1])
    a0 = zeros(Float64,1,lena)
    
    if(lena < 2)
        return a
    end
    
    if(lena == 2)
        if(a[1,1] > a[2,1])
            a0 = a[1,:]
            a[1,:] = a[2,:]
            a[2,:] = a0
        end
        return a
    end
         
    if(lena > 2)
        pivot = a[div(lena,2)]
        i = 1
        j = lena

        while(i <= j)
            
            while(a[i,1] < pivot)
                i = i + 1
            end
            
            while(a[j,1] > pivot)
                j = j - 1
            end

            if(i <= j)
                a0 = a[i,:]
                a[i,:] = a[j,:]
                a[j,:] = a0
                i = i + 1
                j = j - 1
            end
        end

        a[1:i,:] = quicksort2D(a[1:i,:])
        a[i:end,:] = quicksort2D(a[i:end,:]) 

    end
   
    return a
    
end

function run_quicksort()

#Data import
pldata = readdlm("mystery_planet.txt")
time_data = pldata[:,1]
RV_data = pldata[:,2]
err_data = pldata[:,3]

#Sort by phase, given period
phase = Array(Real,length(time_data))
period = 10.
phase = mod(time_data,period)

newsorted1 = [phase RV_data]
data = quicksort2D(newsorted1)

println(data)

end

run_quicksort()