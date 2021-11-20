#! /bin/csh
#
awk '{ 	if( l1 != $1 ){ \
		if( cnt != 0 ){ \
			avg = sum / cnt \
			sd = 0 \
			for( i = 1; i <= cnt; i++ ) \
				sd += ( x[i] - avg ) * ( x[i] - avg ) \
			if( cnt > 1 ) \
				sd = sqrt( sd / ( cnt - 1 ) ) \
			else \
				sd = 0 \
			printf( "%-10s %4d %8.3f %8.3f %8.3f %8.3f\n", l1, cnt, avg, sd, min, max ) \
		} \
		cnt = 0 \
		sum = 0 \
		max = min = $2 \
	} \
	if( $2 > max ) \
		max = $2 \
	if( $2 < min ) \
		min = $2 \
	l1 = $1 \
	cnt++ \
	sum += $2 \
	x[ cnt ] = $2 \
     } \
     END { \
	if( cnt > 0) \
		avg = sum / cnt \
	else \
		avg = 0 \
	sd = 0 \
	for( i = 1; i <= cnt; i++ ) \
		sd += ( x[i] - avg ) * ( x[i] - avg ) \
	if( cnt > 1 ) \
		sd = sqrt( sd / ( cnt - 1 ) ) \
	else \
		sd = 0 \
	if( cnt > 0) \
		printf( "%-10s %4d %8.3f %8.3f %8.3f %8.3f\n", l1, cnt, avg, sd, min, max ) } ' $1
