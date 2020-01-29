classdef Lshape
    
    properties
        b_l(1,2)
        h_l(1,2)
        
    end
    
    methods
        
        function A=a_l(obj)
            A=(obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2));
        end
        
        function Iy=iy(obj)
             Z_bar=((obj.b_l(1).*obj.h_l(1)).*(obj.h_l(2)+obj.h_l(1)./2)+(obj.b_l(2).*obj.h_l(2)).*obj.h_l(2)./2)./((obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2)));
             Iy=(((obj.b_l(1).*obj.h_l(1))./12.*obj.h_l(1).^2)+((obj.b_l(1).*obj.h_l(1)).*(obj.h_l(2)+obj.h_l(1)./2-Z_bar).^2))+(((obj.b_l(2).*obj.h_l(2))./12*obj.h_l(2).^2)+((obj.b_l(2).*obj.h_l(2)).*(obj.h_l(2)./2-Z_bar).^2));
        end
        function Iz=iz(obj)
            Y_bar=((obj.b_l(1).*obj.h_l(1)).*(obj.b_l(1)./2)+(obj.b_l(2).*obj.h_l(2)).*obj.b_l(2)./2)./((obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2)));
            Iz=(((obj.b_l(1).*obj.h_l(1))./12.*obj.b_l(1).^2)+((obj.b_l(1).*obj.h_l(1)).*(obj.b_l(1)./2-Y_bar).^2))+(((obj.b_l(2).*obj.h_l(2))./12.*obj.b_l(2).^2)+((obj.b_l(2).*obj.h_l(2)).*(obj.b_l(2)./2-Y_bar).^2));
        end
        function IYZ=Iyz(obj)
%             Y_bar=((obj.b_l(1).*obj.h_l(1)).*(obj.b_l(1)./2)+(obj.b_l(2).*obj.h_l(2)).*obj.b_l(2)./2)./((obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2)));
%             Z_bar=((obj.b_l(1).*obj.h_l(1)).*(obj.h_l(2)+obj.h_l(1)./2)+(obj.b_l(2).*obj.h_l(2)).*obj.h_l(2)./2)./((obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2)));
            A=(obj.b_l(1).*obj.h_l(1))+(obj.b_l(2).*obj.h_l(2));
            IYZ=((obj.b_l(1).^2./4).*(obj.b_l(2).^2+(obj.h_l(1)+obj.h_l(2)).^2-obj.b_l(1).^2))-(A.*(1./A.*(obj.b_l(1)./2.*(obj.b_l(2).^2+(obj.b_l(1).*(obj.h_l(1)+obj.h_l(2)))-obj.b_l(1).^2))).*(1./A.*(obj.b_l(1)./2.*((obj.h_l(1)+obj.h_l(2)).^2+(obj.b_l(1).*(obj.b_l(2)))-obj.b_l(1).^2))));
            
        end
    end
end
