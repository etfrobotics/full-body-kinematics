function [] = plotJoints(H)
n = length(H);
quiver3(reshape(H(1,4,:),[n,1]), ...
    reshape(H(2,4,:),[n,1]), ...
    reshape(H(3,4,:),[n,1]), ...
    reshape(H(1,1,:),[n,1]), ...
    reshape(H(2,1,:),[n,1]), ...
    reshape(H(3,1,:),[n,1]), ...
    0.2,'r')
hold on
quiver3(reshape(H(1,4,:),[n,1]), ...
    reshape(H(2,4,:),[n,1]), ...
    reshape(H(3,4,:),[n,1]), ...
    reshape(H(1,2,:),[n,1]), ...
    reshape(H(2,2,:),[n,1]), ...
    reshape(H(3,2,:),[n,1]), ...
    0.2,'g')
quiver3(reshape(H(1,4,:),[n,1]), ...
    reshape(H(2,4,:),[n,1]), ...
    reshape(H(3,4,:),[n,1]), ...
    reshape(H(1,3,:),[n,1]), ...
    reshape(H(2,3,:),[n,1]), ...
    reshape(H(3,3,:),[n,1]), ...
    0.2,'b')
plot3(reshape(H(1,4,1:5),[5,1]), ...
    reshape(H(2,4,1:5),[5,1]), ...
    reshape(H(3,4,1:5),[5,1]),color="#0072BD")
plot3(reshape(H(1,4,6:13),[8,1]), ...
    reshape(H(2,4,6:13),[8,1]), ...
    reshape(H(3,4,6:13),[8,1]),color="#0072BD")
plot3(reshape(H(1,4,14:21),[8,1]), ...
    reshape(H(2,4,14:21),[8,1]), ...
    reshape(H(3,4,14:21),[8,1]),color="#0072BD")
plot3(reshape(H(1,4,22:30),[9,1]), ...
    reshape(H(2,4,22:30),[9,1]), ...
    reshape(H(3,4,22:30),[9,1]),color="#0072BD")
plot3(reshape(H(1,4,31:39),[9,1]), ...
    reshape(H(2,4,31:39),[9,1]), ...
    reshape(H(3,4,31:39),[9,1]),color="#0072BD")
plot3(reshape(H(1,4,[22,31]),[2,1]), ...
    reshape(H(2,4,[22,31]),[2,1]), ...
    reshape(H(3,4,[22,31]),[2,1]),color="#0072BD")
plot3(reshape(H(1,4,[6,14]),[2,1]), ...
    reshape(H(2,4,[6,14]),[2,1]), ...
    reshape(H(3,4,[6,14]),[2,1]),color="#0072BD")

xlabel('x')
ylabel('y')
zlabel('z')
axis equal
view(120,20)
axis([-1 1 -1 1 -1 1])

end