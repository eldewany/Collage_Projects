classdef Cshape
    
    properties
        b_c(1,2)
        h_c(1,2)
        
    end
    
    methods
        
        function A=a_c(obj)
            A=(obj.b_c(1).*obj.h_c(1))+(2.*(obj.b_c(2).*obj.h_c(2)));
        end
        
        function Iy=iy(obj)
             Iy=((obj.b_c(1).*obj.h_c(1))./12.*obj.h_c(1).^2)+(2.*((obj.b_c(2).*obj.h_c(2)).*((obj.h_c(1)+obj.h_c(2))./2).^2+((obj.b_c(2).*obj.h_c(2))./12*obj.h_c(2).^2)));
        end
        function Iz=iz(obj)
            Y_bar=((obj.b_c(1).*obj.h_c(1)).*obj.b_c(1)./2+(obj.b_c(2).*obj.h_c(2)).*obj.b_c(2))./((obj.b_c(1).*obj.h_c(1))+2.*(obj.b_c(2).*obj.h_c(2)));
            Iz=((obj.b_c(1).*obj.h_c(1)).*(obj.b_c(1)./2-Y_bar).^2+((obj.b_c(1).*obj.h_c(1))./12.*obj.b_c(1).^2))+(2.*((obj.b_c(2).*obj.h_c(2)).*(obj.b_c(2)./2-Y_bar).^2+((obj.b_c(2).*obj.h_c(2))./12.*obj.b_c(2).^2)));
        end
    end
end
