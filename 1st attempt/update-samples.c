int update_samples() {

    // auxiliar variables
    float minDist = FLT_MAX, dist1 = 0.0f, dist2 = 0.0f;
    int minK = INT_MAX, changes_flag = 0, twoK = (K + 2 - 1) / 2;
    struct spoint p;

    // for each of the samples
    for (int i = 0; i < N; i++) {

        // get current point
        p = RANDOM_SAMPLE[i];

        // default values for minimum calculation
        
        minDist = INT_MAX;
        minK = p.k;

        // for each of the clusters
        for (int j = 0; j < K; j+=2) {

            // calculate the euclidian distance
            dist1 = euclidianDistance(p, (CLUSTERS[j]).centroid);
            dist2 = (j+1 < K) ? (euclidianDistance(p, (CLUSTERS[j+1]).centroid)) : INT_MAX;
            // if a new minimum is found

            if (dist1 <= dist2) {
                if (dist1 < minDist){
                    minDist = dist1;
                    minK = j;
                }
            }
            else {
                if (dist2 < minDist) {
                    minDist = dist2;
                    minK = j+1;
                }
            }
        }


        if (minK != (p.k)) {
            
            // assuming k will always be >= 0
            // update previous cluster
            if (p.k != -1)
                (CLUSTERS[p.k]).dimension--;
            // update newer cluster
            (CLUSTERS[minK]).dimension++;
            // update point's cluster
            (RANDOM_SAMPLE[i]).k = minK;
            changes_flag = 1;
        }
    }

    return changes_flag;
}