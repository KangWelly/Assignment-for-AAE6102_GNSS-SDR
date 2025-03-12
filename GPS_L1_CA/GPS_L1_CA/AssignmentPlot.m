%% for multi-correaltor output
data=trackResults(1).I_multi{200};
if data(6)<0
    data=-data
end
plot(-0.5:0.1:0.5,data,'r-');hold on;
scatter(-0.5:0.1:0.5,data,'bo');
xlabel('code delay');
ylabel('ACF');
title('ACF of Multi-correlator');

%% for Weighted Least Square for positioning (elevation based)
open_gt=[22.328444770087565,114.1713630049711,3];
geobasemap satellite;
error=[];
for i=1:size(navSolutions.latitude,2)
    geoplot(navSolutions.latitude(i),navSolutions.longitude(i),'r*', 'MarkerSize', 10);hold on;
end
  geoplot(open_gt(1),open_gt(2),'o','MarkerFaceColor','y', 'MarkerSize', 10,'MarkerEdgeColor','y');hold on;

%% WLS for velocity
v=[];
for i=2:size(navSolutions.X,2)
    v=[v;norm([navSolutions.X(i),navSolutions.Y(i),navSolutions.Z(i)]-[navSolutions.X(i-1),navSolutions.Y(i-1),navSolutions.Z(i-1)])];
end
plot(1:size(navSolutions.X,2)-1,v);
xlabel('epoch(s)');
ylabel('Velocity(m/s)');


%% for Kalman Filter
